# 04_train_models.R
train_models <- function(config, data, filter_probes_dir, out_dir) {
  
  lgr::lgr$info("Starting model training...")

  train_models_dir <- file.path(out_dir, "04_train_models")
  if (!dir.exists(train_models_dir)) dir.create(train_models_dir, recursive = TRUE)

  pheno_vec <- data$pheno[[config$process_data_params$pheno$pheno_name]]
  methyl_mat <- data$methyl

  filtered_files <- list.files(filter_probes_dir, pattern = "\\.csv$", full.names = TRUE)
  filtered_files <- filtered_files[!grepl("_performance_curve\\.csv$", filtered_files)]
  train_params <- config$train_models_params
  
  models_to_train <- names(train_params$models)
  
  n_cores <- config$main_pipeline_params$n_cores
  cl <- NULL
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
  }

  if (is.null(train_params$tuning_metric)) {
    stop("tuning_metric must be specified in train_models_params")
  }
  tuning_metric <- train_params$tuning_metric
  
  lgr::lgr$info(paste("Tuning models based on metric:", tuning_metric))
  
  # Loop through each probe list
  for (filtered_file in filtered_files) {
    ranking_method <- tools::file_path_sans_ext(basename(filtered_file))
    
    probe_df <- utils::read.csv(filtered_file)
    n_probes <- nrow(probe_df)
    probe_names <- probe_df$probe
    
    lgr::lgr$info(paste0("Training models for probe set: ", ranking_method, " (", n_probes, " probes)"))
    
    method_model_dir <- file.path(train_models_dir, ranking_method)
    if (!dir.exists(method_model_dir)) dir.create(method_model_dir)
    
    train_df <- data.frame(pheno = pheno_vec, t(methyl_mat[probe_names, , drop = FALSE]))
    
    for (model_name in models_to_train) {
      lgr::lgr$info(paste0("Training model: ", model_name))
      model_params <- train_params$models[[model_name]]
      
      set.seed(config$main_pipeline_params$seed)

      if (!is.null(model_params$tuneLength)) {
        model_info <- caret::getModelInfo(model_name)[[model_name]]
        tune_grid_for_len <- model_info$grid(x = train_df[, -1], y = train_df$pheno, len = model_params$tuneLength)
        n_combinations <- nrow(tune_grid_for_len)
      } else {
        stop(paste("For model '", model_name, "', a value for tuneLength is required in the config.", sep=""))
      }
      
      lgr::lgr$info(paste0("Determined ", n_combinations, " tuning combinations for ", model_name))

      num_resamples <- model_params$folds * model_params$repeats
      seeds <- vector(mode = "list", length = num_resamples + 1)
      for(i in 1:num_resamples) seeds[[i]] <- sample.int(100000, n_combinations)
      seeds[[num_resamples + 1]] <- sample.int(100000, 1)
      
      train_control <- caret::trainControl(
        method = "repeatedcv",
        number = model_params$folds,
        repeats = model_params$repeats,
        allowParallel = TRUE, 
        returnData = FALSE,
        savePredictions = "none",
        trim = TRUE,
        seeds = seeds,
        summaryFunction = r_summary
      )

      # Need to turn off parallel processing here for xgboost to work
      if (model_name %in% c("xgbTree", "xgbLinear")) {
        train_control$allowParallel <- FALSE
      }

      train_args <- list(form = pheno ~ ., data = train_df, method = model_name, trControl = train_control,
                         metric = tuning_metric, verbose = FALSE)
      
      train_args$tuneLength <- model_params$tuneLength

      # Setting this for xgboost prevents spam output about a deprecated argument name that is used by caret
      if (model_name %in% c("xgbTree", "xgbLinear")) {
        train_args$verbosity <- 0
      }
      
      model <- do.call(caret::train, train_args)

      best_tune_params <- model$bestTune
      param_values <- sapply(names(best_tune_params), function(p_name) {
        val <- best_tune_params[[p_name]]; if (is.numeric(val)) round(val, 4) else as.character(val)
      })
      params_str <- paste(names(best_tune_params), param_values, sep = " = ", collapse = ", ")
      lgr::lgr$info(paste("Best tune for", model_name, ":", params_str))

      attr(model$terms, ".Environment") <- NULL # Makes the model object smaller, saving space
      fold_id <- as.integer(gsub("fold_", "", basename(out_dir)))
      
      output_object <- list(
        model = model,
        metadata = list(
          ranking_method = ranking_method,
          training_method = model_name,
          n_probes = n_probes,
          outer_fold = fold_id,
          config = config
        )
      )

      if (nrow(model$results) > 1) {
        pdf(file.path(method_model_dir, paste0("f_", fold_id, "_r_", ranking_method, "_m_", model_name, "_tune_plot.pdf")),
            width = 8, height = 6)
        print(plot(model))
        dev.off()
      } else {
        lgr::lgr$info(paste("Skipping tuning plot for", model_name, "as it has only one tuning parameter value."))
      }
      
      filename <- paste0("f_", fold_id, "_r_", ranking_method, "_m_", model_name, ".rds")
      saveRDS(output_object, file.path(method_model_dir, filename))
    }
  }
    
  if (!is.null(cl)) {
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
  }
  
  lgr::lgr$info("Model training complete.")
  return(train_models_dir)
}