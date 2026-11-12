# 07_ensemble_models.R
ensemble_models <- function(config, data_raw, out_dir) {
  lgr::lgr$info("Starting global ensemble creation...")
  set.seed(config$main_pipeline_params$seed)
  
  ensemble_params <- config$ensemble_models_params
  
  ensemble_output_dir <- file.path(out_dir, "ensemble")
  if (!dir.exists(ensemble_output_dir)) dir.create(ensemble_output_dir)
  
  # Get candidate models
  lgr::lgr$info("Loading and filtering candidate base models...")
  perf_file <- file.path(out_dir, "performance", "aggregated_test_results.csv")
  if (!file.exists(perf_file)) stop("Aggregated results file not found.")
  
  all_perf_data <- utils::read.csv(perf_file)
  
  # Filter candidate models by performance threshold
  if (is.null(ensemble_params$selection_metric)) stop("selection_metric missing in config.")
  if (is.null(ensemble_params$min_perf)) stop("min_perf missing in ensemble_models_params.")
  
  sort_metric <- ensemble_params$selection_metric
  min_perf <- as.numeric(ensemble_params$min_perf)
  maximise_sort <- config$metric_maximise_lookup[[sort_metric]]
  
  if (is.null(maximise_sort)) stop(paste("Maximise lookup not found for", sort_metric))
  
  lgr::lgr$info(paste0("Filtering models: ", sort_metric, " must be ", 
                       if (maximise_sort) ">=" else "<=", " ", min_perf))
  
  # Filter models based on threshold
  if (maximise_sort) {
    candidate_pool <- all_perf_data %>% dplyr::filter(.data[[sort_metric]] >= min_perf)
  } else {
    candidate_pool <- all_perf_data %>% dplyr::filter(.data[[sort_metric]] <= min_perf)
  }
  
  # Sort
  candidate_pool <- candidate_pool %>% 
    dplyr::arrange(if (maximise_sort) dplyr::desc(.data[[sort_metric]]) else .data[[sort_metric]])
  
  n_candidates <- nrow(candidate_pool)
  lgr::lgr$info(paste0("Found ", n_candidates, " models that meet the minimum performance threshold."))
  
  if (n_candidates == 0) stop("No models met the minimum performance threshold.")
  
  candidate_paths <- candidate_pool$path
  candidate_names <- tools::file_path_sans_ext(basename(candidate_paths))
  
  # Prepare meta-model data
  lgr::lgr$info("Generating predictions for candidate base models...")
  
  # Load models
  if(any(is.na(candidate_paths) | candidate_paths == "")){ stop("Invalid model paths.") }
  pool_models <- lapply(candidate_paths, readRDS)
  names(pool_models) <- candidate_names
  
  # Process full dataset for meta-model
  metamodel_process_dir <- file.path(ensemble_output_dir, "metamodel_processing")
  process_data_result <- process_data(config = config, data = data_raw, out_dir = metamodel_process_dir, process_methyl = TRUE, process_pheno = TRUE)
  processed_data <- readRDS(file.path(process_data_result$dir, "processed_data.rds"))
  
  # Generate predictions
  preds_list <- lapply(pool_models, function(bm) {
    model <- bm$model
    probes_needed <- model$coefnames
    samples_available <- colnames(processed_data$methyl)
    
    if (!all(probes_needed %in% rownames(processed_data$methyl))) return(NULL)
    
    methyl_subset_t <- t(processed_data$methyl[probes_needed, samples_available, drop = FALSE])
    pred_df <- as.data.frame(methyl_subset_t)
    return(predict(model, newdata = pred_df))
  })
  
  # Clean nulls
  valid_indices <- !sapply(preds_list, is.null)
  preds_list <- preds_list[valid_indices]
  pool_models <- pool_models[valid_indices]
  candidate_names <- candidate_names[valid_indices]
  n_candidates <- length(candidate_names)
  
  X_full <- as.data.frame(do.call(cbind, preds_list))
  y_full <- processed_data$pheno[[config$process_data_params$pheno$pheno_name]]
  
  # Search for optimal combination of base models
  lgr::lgr$info("Using PRESS to identify optimal selection of base models...")
  
  max_k <- if (!is.null(ensemble_params$max_models)) as.integer(ensemble_params$max_models) else 3 # Define range of ensemble sizes (1 to max_models)
  
  # Ensure max_k is not larger than the number of available candidate models
  if (max_k > n_candidates) {
    lgr::lgr$warn(paste0("max_models (", max_k, ") exceeds available candidates (", n_candidates, "). Adjusting maximum base model count."))
    max_k <- n_candidates
  }
  
  lgr::lgr$info(paste0("Searching combinations between 1 and ", max_k, "..."))
  
  TSS <- sum((y_full - mean(y_full))^2) # Get total sum of squares
  
  # Get PRESS statistic and predicted R squared
  calc_press_r2 <- function(X_subset, y, TSS) {
    dat <- as.data.frame(X_subset)
    dat$target <- y
    
    fit <- try(stats::lm(target ~ ., data = dat), silent = TRUE)
    if (inherits(fit, "try-error")) return(-Inf)
    
    h <- stats::hatvalues(fit)
    if (any(h >= 1 - 1e-10)) return(-Inf)
    
    press_residuals <- fit$residuals / (1 - h)
    PRESS <- sum(press_residuals^2)
    
    # Predicted R-squared (similar to 1 - RSS/TSS)
    pred_r2 <- 1 - (PRESS / TSS)
    return(pred_r2)
  }
  
  best_r2 <- -Inf
  best_combo_names <- NULL
  
  all_results_list <- list()
  counter <- 1
  
  # Go through ensemble sizes
  for (k in 1:max_k) {
    combos <- combn(candidate_names, k)
    lgr::lgr$info(paste0("Evaluating ", ncol(combos), " combinations of size ", k, "..."))
    
# Go through specific combinations of given size
for (i in 1:ncol(combos)) {
  current_names <- combos[, i]
  X_sub <- X_full[, current_names, drop = FALSE]
  
  if (k == 1) { # For single models, get R squared directly from base model predictions
    preds <- X_sub[[1]]
    SS_res <- sum((y_full - preds)^2)
    current_r2 <- 1 - (SS_res / TSS)
  } else {
    current_r2 <- calc_press_r2(X_sub, y_full, TSS)
  }
  
  all_results_list[[counter]] <- list(
        combination_id = counter,
        n_models = k,
        models = paste(current_names, collapse = " | "),
        loocv_pred_r2 = current_r2
      )
      
      # Identify best selection so far
      if (!is.na(current_r2) && current_r2 > best_r2) {
        best_r2 <- current_r2
        best_combo_names <- current_names
      }
      
      counter <- counter + 1
    }
  }
  
  # Save results
  all_combos_df <- dplyr::bind_rows(all_results_list) %>%
    dplyr::arrange(dplyr::desc(loocv_pred_r2))
  
  combos_csv_path <- file.path(ensemble_output_dir, "ensemble_combinations_performance.csv")
  utils::write.csv(all_combos_df, combos_csv_path, row.names = FALSE)
  lgr::lgr$info(paste("Saved performance of all", nrow(all_combos_df), "combinations to:", combos_csv_path))
  
  lgr::lgr$info(paste0("Best ensemble selected, with predicted R squared of ", round(best_r2, 4)))
  lgr::lgr$info(paste("Models (", length(best_combo_names), "): ", paste(best_combo_names, collapse = ", ")))
  
  # Finalise ensemble
  lgr::lgr$info("Generating final meta-model...")
  
  final_base_models <- pool_models[best_combo_names]
  meta_training_data <- X_full[, best_combo_names, drop = FALSE]
  meta_training_data$target <- y_full
  
  meta_model <- stats::lm(target ~ ., data = meta_training_data)
  meta_summary <- summary(meta_model)
  
  global_ensemble <- list(
    meta_model = meta_model,
    base_models = final_base_models,
    metadata = list(
      config = config,
      selected_models = best_combo_names,
      loocv_metric = "Predicted R squared",
      loocv_score = best_r2
    )
  )
  
  pheno_name <- config$process_data_params$pheno$pheno_name
  
  # Make model object smaller
  global_ensemble$meta_model$call <- quote(dummy_call())
  global_ensemble$meta_model$residuals <- NULL
  global_ensemble$meta_model$fitted.values <- NULL
  global_ensemble$meta_model$model <- NULL
  attr(global_ensemble$meta_model$terms, ".Environment") <- baseenv()
  for (i in seq_along(global_ensemble$base_models)) {
    attr(global_ensemble$base_models[[i]]$model$terms, ".Environment") <- baseenv()
  }
  for (i in seq_along(global_ensemble$base_models)) {
    global_ensemble$base_models[[i]]$model$call <- quote(dummy_call())
  }
  
  output_filename <- paste0("ens_", pheno_name, ".rds")
  output_path <- file.path(ensemble_output_dir, output_filename)
  saveRDS(global_ensemble, output_path)
  
  lgr::lgr$info(paste("Saved ensemble to:", output_path))
  lgr::lgr$info("Ensemble creation complete.")
  
  return(output_path)
}