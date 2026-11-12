# 03_select_probes.R
select_probes <- function(config, data, rank_probes_dir, out_dir) {
  
  lgr::lgr$info("Starting probe selection...")
  
  select_probes_dir <- file.path(out_dir, "03_select_probes")
  if (!dir.exists(select_probes_dir)) dir.create(select_probes_dir, recursive = TRUE)
  
  train_pheno_vec <- data$pheno[[config$process_data_params$pheno$pheno_name]]
  train_methyl_mat <- data$methyl
  validation_pheno_vec <- train_pheno_vec
  validation_methyl_mat <- train_methyl_mat
  
  ranking_files <- list.files(rank_probes_dir, pattern = "ranking_.*\\.csv", full.names = TRUE)
  select_params <- config$select_probes_params
  
  n_cores <- config$main_pipeline_params$n_cores
  cl <- NULL
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
  }
  
  # Loop through each ranking method
  for (rank_file in ranking_files) {
    method_name <- gsub("ranking_|\\.csv", "", basename(rank_file))
    
    if (is.null(select_params$probe_counts[[method_name]])) {
      lgr::lgr$info(paste0("Skipping probe selection for '", method_name, "' as it lacks a config in 'select_probes_params$probe_counts'."))
      next
    }
    
    lgr::lgr$info(paste("Selecting probes for ranking method:", method_name))
    ranked_probes_df_full <- utils::read.csv(rank_file)
    final_probe_counts_for_method <- c()
    
    # Performance-based selection
      lgr::lgr$info("Starting performance-based selection...")
      perf_params <- select_params
      
      max_probes_config_val <- perf_params$max_probes[[method_name]]
      if (is.null(max_probes_config_val)) {
        lgr::lgr$warn(paste0("max_probes not specified for '", method_name, "'. Using all ", nrow(ranked_probes_df_full), " probes as the upper limit."))
        max_probes_initial <- nrow(ranked_probes_df_full)
      } else {
        max_probes_initial <- max_probes_config_val
      }

      ranked_probes_for_selection <- ranked_probes_df_full %>% dplyr::slice(1:min(max_probes_initial, n()))
      max_probes_final <- nrow(ranked_probes_for_selection)
      n_probes_seq <- round(seq(from = perf_params$min, to = max_probes_final, length.out = perf_params$n_counts))
      n_probes_seq <- n_probes_seq[n_probes_seq <= max_probes_final]
      
      if (length(n_probes_seq) == 0) {
        lgr::lgr$warn("Probe sequence is empty. Skipping performance selection.")
        next
      }
      
      lgr::lgr$info(paste("Initial probe counts for performance curve:", paste(n_probes_seq, collapse = ", ")))
      
      if (!is.null(select_params$corr_filter_thresh)) {
        lgr::lgr$info("Pre-calculating final probe counts after correlation filtering...")
        filtered_counts_df <- foreach::foreach(n = n_probes_seq, .combine = 'rbind', .packages = 'caret', .export = "r_summary") %dopar% {
          if (n < 2) {
            return(data.frame(n_initial = n, n_final = n))
          }
          probes_to_consider_df <- ranked_probes_for_selection[1:n, ]
          sub_methyl_mat <- train_methyl_mat[probes_to_consider_df$probe, , drop = FALSE]
          cor_matrix <- cor(t(sub_methyl_mat))
          probes_to_remove_idx <- caret::findCorrelation(cor_matrix, cutoff = select_params$corr_filter_thresh, exact = FALSE)
          final_count <- n - length(probes_to_remove_idx)
          data.frame(n_initial = n, n_final = final_count)
        }
        count_summary <- paste(filtered_counts_df$n_initial, "->", filtered_counts_df$n_final, collapse = ", ")
        lgr::lgr$info(paste("Probe counts to be tested (initial -> post-filtering):", count_summary))
      }
      
      set.seed(config$main_pipeline_params$seed)
      tuning_metric <- perf_params$tuning_metric
      
      selection_metric <- if (!is.null(perf_params$selection_metric)) perf_params$selection_metric else tuning_metric
      lgr::lgr$info(paste0("Using '", tuning_metric, "' for hyperparameter tuning."))
      lgr::lgr$info(paste0("Using '", selection_metric, "' for performance curve and probe count selection."))
    
      maximise <- config$metric_maximise_lookup[[selection_metric]]
      model_config <- perf_params$models[[method_name]]
      
      if (is.null(model_config)) {
        model_config <- c(method = "ranger", folds = 5, repeats = 1, tuneLength = 5)
        lgr::lgr$info(paste("Model for '", method_name, "' is NULL, so using default: ranger.", sep = ""))
      }
      
      model_name <- model_config['method']
      model_params <- model_config
      
      lgr::lgr$info(paste0("Generating performance curve for '", method_name, "' using the '", model_name, "' model..."))
      
      train_control <- caret::trainControl(method = "repeatedcv", number = as.numeric(model_params["folds"]), repeats = as.numeric(model_params["repeats"]), summaryFunction = r_summary)
      
      results <- foreach::foreach(n = n_probes_seq, .combine = 'rbind', .packages = c('caret', 'lgr', 'dplyr'), .export = "r_summary") %dorng% {
        top_n_probes_df <- ranked_probes_for_selection[1:n, ]
        final_probes_for_training_df <- top_n_probes_df
        if (!is.null(select_params$corr_filter_thresh) && nrow(top_n_probes_df) >= 2) {
          sub_methyl_mat <- train_methyl_mat[top_n_probes_df$probe, , drop = FALSE]
          cor_matrix <- cor(t(sub_methyl_mat))
          probes_to_remove_idx <- caret::findCorrelation(cor_matrix, cutoff = select_params$corr_filter_thresh, exact = FALSE)
          if (length(probes_to_remove_idx) > 0) {
            final_probes_for_training_df <- top_n_probes_df[-probes_to_remove_idx, ]
          }
        }
        final_probes <- final_probes_for_training_df$probe
        train_df <- data.frame(pheno = train_pheno_vec, t(train_methyl_mat[final_probes, , drop = FALSE]))
        validation_df <- data.frame(pheno = validation_pheno_vec, t(validation_methyl_mat[final_probes, , drop = FALSE]))
        
        train_args <- list(form = pheno ~ ., data = train_df, method = model_name, trControl = train_control, metric = tuning_metric, verbose = FALSE)
        train_args$tuneLength <- as.numeric(model_params["tuneLength"])
        
        model <- do.call(caret::train, train_args)
        predictions <- predict(model, newdata = validation_df)
        eval_data <- data.frame(pred = predictions, obs = validation_df$pheno)
        
        perf_orig <- r_summary(eval_data)[selection_metric]
        data.frame(n_probes = n, n_probes_final = length(final_probes), performance = perf_orig)
      }
      
      utils::write.csv(results, file.path(select_probes_dir, paste0(method_name, "_performance_curve.csv")), row.names = FALSE)
      
      selection_configs <- perf_params$probe_counts[[method_name]]
      optimal_probes_df <- data.frame(n_probes = integer())
      
      for (config_item in selection_configs) {
        
      # Get values for the given range and the specified number of probe sets
      range_bounds <- as.numeric(config_item[["range"]])
      counts <- as.integer(config_item["counts"])
      
      # Define probe count bounds from the config
      lower_bound <- range_bounds[1]
      upper_bound <- range_bounds[2]
        
        # Filter the results df
        filtered_results <- results %>%
          dplyr::filter(n_probes >= lower_bound & n_probes <= upper_bound)
        
        if (nrow(filtered_results) > 0) {
          sorted_results <- if (maximise) {
            filtered_results %>% dplyr::arrange(desc(performance))
          } else {
            filtered_results %>% dplyr::arrange(performance)
          }
          
          selected_probes <- sorted_results %>%
            dplyr::slice(1:min(counts, n())) %>%
            dplyr::pull(n_probes)
          
          if (length(selected_probes) > 0) {
            optimal_probes_df <- rbind(optimal_probes_df, data.frame(n_probes = selected_probes))
          }
        }
      }
      
      optimal_probes_df <- optimal_probes_df %>% dplyr::distinct()
      
      final_probe_counts_for_method <- c(final_probe_counts_for_method, unique(optimal_probes_df$n_probes))
      
      p <- ggplot2::ggplot(results, ggplot2::aes(x = n_probes, y = performance)) +
        ggplot2::geom_line(linewidth = 1, color = "blue") +
        ggplot2::labs(
          title = paste("Performance curve for", method_name), 
          subtitle = paste("Model:", model_name),
          x = "Number of probes (before filtering)", 
          y = selection_metric
        ) +
        ggplot2::theme_minimal()
      if (nrow(optimal_probes_df) > 0) {
        p <- p + ggplot2::geom_vline(data = optimal_probes_df, ggplot2::aes(xintercept = n_probes), linetype = "dashed")
      }
      ggplot2::ggsave(file.path(select_probes_dir, paste0(method_name, "_perf_plot.pdf")), p, width = 8, height = 6)
    
    # Save final probe lists
    final_probe_counts <- unique(final_probe_counts_for_method)
    lgr::lgr$info(paste0("Final probe counts for ", method_name, " (before filtering): ", paste(final_probe_counts, collapse = ", ")))
    for (n_probes in final_probe_counts) {
      if (n_probes > 0 && n_probes <= nrow(ranked_probes_df_full)) {
        probes_to_consider_df <- ranked_probes_df_full[1:n_probes, ]
        final_probes_df <- probes_to_consider_df
        if (!is.null(select_params$corr_filter_thresh)) {
          lgr::lgr$info(paste0("Applying correlation filter to probe set '", method_name, "' of size ", n_probes))
          probes_exist <- probes_to_consider_df$probe %in% rownames(train_methyl_mat)
          if (!all(probes_exist)) {
            lgr::lgr$warn(paste0(sum(!probes_exist), " probes not found in methylation data for this set. Filtering on available probes."))
            probes_to_consider_df <- probes_to_consider_df[probes_exist, ]
          }
          if (nrow(probes_to_consider_df) < 2) {
            lgr::lgr$info("Skipping correlation filtering as there are fewer than two probes.")
          } else {
            cor_matrix <- cor(t(train_methyl_mat[probes_to_consider_df$probe, , drop = FALSE]))
            probes_to_remove_idx <- caret::findCorrelation(cor_matrix, cutoff = select_params$corr_filter_thresh, exact = FALSE)
            if (length(probes_to_remove_idx) > 0) {
              probes_to_remove_names <- rownames(cor_matrix)[probes_to_remove_idx]
              final_probes_df <- probes_to_consider_df %>%
                dplyr::filter(!probe %in% probes_to_remove_names)
              lgr::lgr$info(paste0(length(probes_to_remove_idx), " probes removed, ", nrow(final_probes_df), " retained."))
            } else {
              lgr::lgr$info("No probes were removed by the correlation filter for this set.")
            }
          }
        }
        output_filename <- file.path(select_probes_dir, paste0(method_name, "_", nrow(final_probes_df), ".csv"))
        utils::write.csv(final_probes_df, output_filename, row.names = FALSE)
        lgr::lgr$info(paste0("Saved final probe set (", nrow(final_probes_df), " probes) from initial set of ", n_probes, " to: ", basename(output_filename)))
      }
    }
  }
  
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }
  
  lgr::lgr$info("Probe selection complete.")
  return(select_probes_dir)
}