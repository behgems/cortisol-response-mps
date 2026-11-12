# 04_analyse_shap_importance.R
# May need to force low-level libraries to use 1 thread, to prevent nested parallel thrashing. If it hangs at 100% CPU usage, try uncommenting this:
# Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

library(fastshap)
library(shapviz)
library(ggplot2)
library(dplyr)
library(purrr)
library(patchwork)
library(tibble)
library(tidyr)
library(future)
library(furrr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(caret) 
library(viridis)
library(stringr)
library(readr)

n_top_probes <- 10
n_shap_sim <- 1000
n_cores <- 12 # 5

options(future.globals.maxSize = 8 * 1024^3)

ens <- readRDS("./output/02_generate_mps/ensemble/ens_cortisol_auc.rds")
main_out_dir <- file.path("./output/04_analyse_shap_importance")
champagne_pre_pred_proc <- readRDS("./output/03_generate_predictions/champagne_pre_pred/ens_cortisol_auc_on_champagne_pre_pred_proc.rds")
champagne_post_pred_proc <- readRDS("./output/03_generate_predictions/champagne_post_pred/ens_cortisol_auc_on_champagne_post_pred_proc.rds")
GSE77445_pred_proc <- readRDS("./output/03_generate_predictions/GSE77445_pred/ens_cortisol_auc_on_GSE77445_pred_proc.rds")

dir.create(main_out_dir, showWarnings = FALSE, recursive = TRUE)

pred_wrapper <- function(model, newdata) {
  predict(model, as.data.frame(newdata))
}

run_shap_pipeline <- function(champ_pre_proc, champ_post_proc, gse_proc) {
  
  set.seed(456)

  # Only the probes any base model actually uses
  all_model_probes <- unique(unlist(purrr::map(ens$base_models, ~ .x$model$finalModel$xNames)))
  base_models_slim <- purrr::map(ens$base_models, "model")

  datasets_to_explain <- list(
    GSE77445       = t(gse_proc$methyl)[, all_model_probes, drop = FALSE],
    Champagne_Pre  = t(champ_pre_proc$methyl)[, all_model_probes, drop = FALSE],
    Champagne_Post = t(champ_post_proc$methyl)[, all_model_probes, drop = FALSE]
  )
  
  message("Reconstructing training data for use as background data...")
  
  # Get all sample IDs from the processed GSE77445 data
  all_sample_ids <- colnames(gse_proc$methyl)
  rm(champ_pre_proc, champ_post_proc, gse_proc)
  gc()
  
  # Get folds from model object
  folds <- ens$metadata$config$outer_loop_params$folds
  
  # Identify which folds were used to train the base models that are included in the ensemble
  model_fold_keys <- unique(stringr::str_extract(names(ens$base_models), "f_\\d+"))
  fold_numbers <- as.numeric(gsub("f_", "", model_fold_keys))
  
  # Create the background data list by subsetting the full training data set
  background_training_data <- purrr::map(fold_numbers, function(i) {
    fold_name_in_list <- paste0("Fold", i)
    test_indices <- folds[[fold_name_in_list]]
    test_ids <- all_sample_ids[test_indices]
    
    train_ids <- setdiff(all_sample_ids, test_ids)
    
    # Subset transposed matrix (samples are rows)
    datasets_to_explain$GSE77445[train_ids, ]
  }) %>%
    purrr::set_names(model_fold_keys)
  
  # Get SHAP values for each base model
  message("Calculating SHAP values...")
  
  plan(multisession, workers = n_cores)
  
  base_model_names <- names(ens$base_models)
  meta_weights <- coef(ens$meta_model)
  
  computation_grid <- crossing(
    dataset_name = names(datasets_to_explain),
    model_name = base_model_names
  )
  
  shap_results_flat <- future_pmap(computation_grid, function(dataset_name, model_name) {
    message(paste0("Processing: ", model_name, " on data set ", dataset_name))
    
    current_data <- datasets_to_explain[[dataset_name]]
    base_model <- base_models_slim[[model_name]]
    
    # Extract the fold key (e.g., "f_1") from the model name to find the correct background data set
    fold_key <- stringr::str_extract(model_name, "f_\\d+")
    
    if (is.na(fold_key) || !fold_key %in% names(background_training_data)) {
        stop(paste("Could not find matching background data for model:", model_name))
    }
    
    background_data <- background_training_data[[fold_key]]
    model_features <- base_model$finalModel$xNames
    
    # Ensure both datasets have the same features in the same order
    data_to_explain_subset <- current_data[, model_features]
    background_data_subset <- background_data[, model_features]

    explain_args <- list(
      object = base_model,
      X = background_data_subset,
      feature_names = model_features,
      pred_wrapper = pred_wrapper,
      newdata = data_to_explain_subset,
      nsim = n_shap_sim,
      adjust = TRUE
    )
    
    shap_matrix <- do.call(fastshap::explain, explain_args)
    return(shap_matrix)
  }, .options = furrr_options(seed = TRUE, packages = c("caret", "stringr")))
  
  # Reconstruct the nested list structure from the parallel results
  all_shap_values <- list()
  for (i in 1:nrow(computation_grid)) {
    dataset <- computation_grid$dataset_name[i]
    model <- computation_grid$model_name[i]
    if (is.null(all_shap_values[[dataset]])) { all_shap_values[[dataset]] <- list() }
    all_shap_values[[dataset]][[model]] <- shap_results_flat[[i]]
  }
  
  message("Aggregating SHAP values using meta-model weights...")
  all_probes <- unique(unlist(map(ens$base_models, ~ .x$model$finalModel$xNames)))
  final_shap_results <- imap(all_shap_values, function(model_shap_list, dataset_name) {
    final_shap_matrix <- matrix(0,
                                nrow = nrow(datasets_to_explain[[dataset_name]]),
                                ncol = length(all_probes),
                                dimnames = list(rownames(datasets_to_explain[[dataset_name]]), all_probes))
    for (model_name in base_model_names) {
      model_weight <- meta_weights[model_name]
      model_shap <- model_shap_list[[model_name]]
      model_probes <- colnames(model_shap)
      final_shap_matrix[, model_probes] <- final_shap_matrix[, model_probes] + (model_shap * model_weight)
    }
    return(final_shap_matrix)
  })
  
  message("Saving raw and aggregated SHAP matrices...")
  saveRDS(all_shap_values, file.path(main_out_dir, paste0("all_shap_values_raw.rds")))
  saveRDS(final_shap_results, file.path(main_out_dir, paste0("final_shap_results_aggregated.rds")))
  
  calculate_shap_metrics <- function(ensemble_shap_matrix, feature_data, dataset_name) {
    
    feature_data_subset <- feature_data[, colnames(ensemble_shap_matrix)]
    sv <- shapviz::shapviz(ensemble_shap_matrix, X = feature_data_subset)
    
    ensemble_imp_vector <- colMeans(abs(sv$S))
    ensemble_global_imp_df <- data.frame(
      probe = names(ensemble_imp_vector),
      mean_abs_shap = unname(ensemble_imp_vector)
    )
    
    message(paste("Getting global interaction score for all probes in data set", dataset_name))
    all_probes_in_sv <- colnames(sv$S) 
    
    global_interaction_strengths <- furrr::future_map_dbl(all_probes_in_sv, function(v_probe) {
      interaction_vec <- shapviz::potential_interactions(sv, v = v_probe) 
      sum(interaction_vec, na.rm = TRUE) 
    }, .options = furrr_options(seed = TRUE))
    
    global_interaction_df <- tibble::tibble(
      probe = all_probes_in_sv,
      global_interaction_strength = global_interaction_strengths
    )
    
    message(paste("Finished calculations for data set", dataset_name))
    
    return(list(
      importance = ensemble_global_imp_df,
      interactions = global_interaction_df
    ))
  }
  
  results_list <- imap(datasets_to_explain, function(feature_data, dataset_name) {
    calculate_shap_metrics(
      ensemble_shap_matrix = final_shap_results[[dataset_name]],
      feature_data = feature_data,
      dataset_name = dataset_name
    )
  })
  
  message("Saving SHAP metrics...")
  saveRDS(results_list, file.path(main_out_dir, paste0("shap_metrics_results_list.rds")))
  
  message("Comparing global importance across data sets using scaled scores...")
  all_global_importance <- purrr::map(results_list, "importance")
  
  scaled_importance <- all_global_importance %>%
    purrr::map(~ .x %>% 
                 dplyr::mutate(scaled_importance = (mean_abs_shap / sum(mean_abs_shap)) * 100) %>% 
                 dplyr::select(probe, scaled_importance))
  
  comparison_df <- scaled_importance %>%
    purrr::reduce(dplyr::full_join, by = "probe")
  
  colnames(comparison_df) <- c("probe", names(all_global_importance))
  
  comparison_df <- comparison_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      mean_scaled_importance = mean(dplyr::c_across(dplyr::where(is.numeric)), na.rm = TRUE),
      sd_scaled_importance = sd(dplyr::c_across(dplyr::where(is.numeric)), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(mean_scaled_importance))
  
  write.csv(comparison_df, file.path(main_out_dir, paste0("shap_global_importance_comparison.csv")), row.names = FALSE)
  
  # Global importance plot
  top_probes_mean <- head(comparison_df, n_top_probes)
  
  p_comp <- ggplot2::ggplot(
    top_probes_mean,
    ggplot2::aes(x = mean_scaled_importance, y = stats::reorder(probe, mean_scaled_importance))
  ) +
    ggplot2::geom_col(fill = "#FCA50A") +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        xmin = mean_scaled_importance - sd_scaled_importance,
        xmax = mean_scaled_importance + sd_scaled_importance
      ),
      height = 0.2,
      orientation = "y"
    ) +
    ggplot2::labs(
      x = "Mean Importance (% Contribution)",
      y = "Probe ID"
    ) +
    ggplot2::theme_bw(base_size = 16)
  
  ggplot2::ggsave(
    file.path(main_out_dir, paste0("shap_comparison_top_plot.pdf")),
    p_comp,
    width = 10,
    height = 8
  )
  
  message("Comparing global interactions across data sets...")
  all_global_interactions <- purrr::map(results_list, "interactions")
  
  # Rename global_interaction_strength col in each df to be specific to the data set
  renamed_global_interactions <- purrr::imap(all_global_interactions, ~{
    .x %>% dplyr::rename(!! .y := global_interaction_strength)
  })
  
  # Join by probe and calculate mean and SD
  global_interaction_comparison_df <- renamed_global_interactions %>%
    purrr::reduce(dplyr::full_join, by = "probe") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      mean_global_interaction = mean(dplyr::c_across(dplyr::where(is.numeric)), na.rm = TRUE),
      sd_global_interaction = sd(dplyr::c_across(dplyr::where(is.numeric)), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(mean_global_interaction))
  
  write.csv(
    global_interaction_comparison_df,
    file.path(main_out_dir, paste0("shap_global_interactions_comparison.csv")),
    row.names = FALSE
  )
  
  plan(sequential)
}

run_shap_pipeline(
  champ_pre_proc = champagne_pre_pred_proc,
  champ_post_proc = champagne_post_pred_proc,
  gse_proc = GSE77445_pred_proc
)