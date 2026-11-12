# 02_rank_probes.R
rank_probes <- function(config, data, out_dir) {
  
  lgr::lgr$info("Starting probe ranking...")
  
  rank_probes_dir <- file.path(out_dir, "02_rank_probes")
  if (!dir.exists(rank_probes_dir)) dir.create(rank_probes_dir, recursive = TRUE)
  
  train_pheno_vec <- data$pheno[[config$process_data_params$pheno$pheno_name]]
  train_methyl_mat <- data$methyl
  
  params <- config$rank_probes_params
  
  prefilter_probes <- function(methyl_matrix, pheno_vector, n_probes_to_keep, method, config, ranking_df = NULL) {
    if (is.null(n_probes_to_keep)) {
      lgr::lgr$info("No pre-filtering requested. Using all probes.")
      return(methyl_matrix)
    }
    
    if (is.null(method)) {
      lgr::lgr$warn("Prefilter method not specified, defaulting to 'correlation'.")
      method <- "correlation"
    }
    
    n_probes_to_filter <- min(n_probes_to_keep, nrow(methyl_matrix))

    if (!is.null(ranking_df)) {
      # Ensure the ranking_df is sorted correctly before selecting top probes
      if (!"rank" %in% colnames(ranking_df) || is.unsorted(ranking_df$rank)) {
        ranking_df <- ranking_df %>% 
          dplyr::arrange(dplyr::desc(raw_importance)) %>%
          dplyr::mutate(rank = dplyr::row_number())
      }
      top_probes <- ranking_df$probe[1:n_probes_to_filter]
      return(methyl_matrix[top_probes, , drop = FALSE])
    }
    
    if (method == "correlation") {
      correlations <- cor(t(methyl_matrix), pheno_vector, use = "pairwise.complete.obs")
      abs_cor <- abs(correlations[,1])
      top_probes_indices <- order(abs_cor, decreasing = TRUE)[1:n_probes_to_filter]
      
    } else if (method == "mutual_information") {
      lgr::lgr$info("Calculating mutual information for all probes...")
      
      n_cores <- config$main_pipeline_params$n_cores
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      
      disc_phenotype <- infotheo::discretize(pheno_vector)
      
      mi_scores <- foreach::foreach(
        i = 1:nrow(methyl_matrix), 
        .combine = "c", 
        .packages = "infotheo"
      ) %dopar% {
        infotheo::mutinformation(
          infotheo::discretize(methyl_matrix[i,]),
          disc_phenotype,
          method = "emp"
        )
      }
      
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
      
      top_probes_indices <- order(mi_scores, decreasing = TRUE)[1:n_probes_to_filter]
      
    } else {
      stop(paste("Invalid prefilter_method specified:", method))
    }
    
    methyl_mat_filtered <- methyl_matrix[top_probes_indices, , drop = FALSE]
    
    return(methyl_mat_filtered)
  }
  
  save_ranking <- function(df, filename, rank_dir, method_name) {
    df_sorted <- df %>% 
      dplyr::arrange(dplyr::desc(raw_importance)) %>%
      dplyr::mutate(rank = dplyr::row_number()) %>%
      dplyr::select(probe, rank, raw_importance)
    
    output_path <- file.path(rank_dir, filename)
    utils::write.csv(df_sorted, output_path, row.names = FALSE)
    lgr::lgr$info(paste("Ranking file saved to:", output_path))
  }
  
  # Run selected ranking methods
  lgr::lgr$info("Ranking probes by absolute correlation...")
  correlations <- cor(t(train_methyl_mat), train_pheno_vec, use = "pairwise.complete.obs")
  corr_df <- data.frame(probe = rownames(train_methyl_mat), raw_importance = abs(correlations[,1]))
  save_ranking(corr_df, "ranking_correlation.csv", rank_probes_dir, "correlation")
  
  lgr::lgr$info("Ranking probes using mutual information...")
    
    n_cores <- config$main_pipeline_params$n_cores
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    lgr::lgr$info("Discretising the phenotype...")
    disc_phenotype <- infotheo::discretize(train_pheno_vec)
    
    probe_indices <- 1:nrow(train_methyl_mat)
    probe_batches <- split(probe_indices, ceiling(seq_along(probe_indices) / 1000))
    
    lgr::lgr$info(paste("Processing", length(probe_indices), "probes in", length(probe_batches), "batch(es)..."))

    mi_scores <- foreach::foreach(
        batch_idx = probe_batches, 
        .combine = "c", 
        .packages = "infotheo"
    ) %dopar% {
    
      batch_data <- train_methyl_mat[batch_idx, , drop = FALSE]
      
      sapply(1:nrow(batch_data), function(i) {
        infotheo::mutinformation(
          infotheo::discretize(batch_data[i,]),
          disc_phenotype,
          method = "emp"
        )
      })
    }

    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    
    mi_df <- data.frame(
        probe = rownames(train_methyl_mat),
        raw_importance = mi_scores
    )
    
    save_ranking(mi_df, "ranking_mutual_information.csv", rank_probes_dir, "mutual_information")

  resolve_prefilter_probes <- function(stab_params, methyl_mat, pheno_vec, config, rank_dir) {
    
    prefilter_params <- stab_params$prefilter_probes
    prefilter_type <- prefilter_params$type
    prefilter_n <- as.numeric(prefilter_params$n)
    prefilter_method <- prefilter_params$method

    lgr::lgr$info(paste("Prefiltering with method '", prefilter_method, "' (n = ", prefilter_n, ")", sep=""))

    ranking_file_path <- file.path(rank_dir, paste0("ranking_", prefilter_method, ".csv"))
    if (!file.exists(ranking_file_path)) {
        stop(paste("Could not find required ranking file for pre-filtering:", ranking_file_path))
    }
    lgr::lgr$info(paste("Loading existing ranking from:", basename(ranking_file_path)))
    prelim_ranking_df <- utils::read.csv(ranking_file_path)

    if (!"rank" %in% colnames(prelim_ranking_df)) {
      prelim_ranking_df <- prelim_ranking_df %>%
        dplyr::arrange(dplyr::desc(raw_importance)) %>%
        dplyr::mutate(rank = dplyr::row_number())
    }

    n_probes_final <- NULL
    knee_point_data <- NULL

    n_probes_final <- floor(prefilter_n)
    lgr::lgr$info(paste("Number of probes to use for pre-filtering:", n_probes_final))

    return(list(
        n_probes = n_probes_final,
        prelim_ranking_df = prelim_ranking_df
    ))
  }
  
if (params$lasso_stability_selection$use) {
    stab_params <- params$lasso_stability_selection
    lgr::lgr$info("Ranking probes using glmnet stability selection...")

    prefilter_resolution <- resolve_prefilter_probes(stab_params, train_methyl_mat, train_pheno_vec, config, rank_probes_dir)
    n_probes_to_keep <- prefilter_resolution$n_probes

    methyl_mat_filtered <- prefilter_probes(
      train_methyl_mat, train_pheno_vec,
      n_probes_to_keep, stab_params$prefilter_probes$method, config,
      ranking_df = prefilter_resolution$prelim_ranking_df
    )

    x_predictors <- t(methyl_mat_filtered)
    y_outcome <- train_pheno_vec

    if (is.null(stab_params$tuning$metric)) {
      stop("Metric must be specified in config.")
    }
    lambda_metric <- stab_params$tuning$metric

    train_control_tuning <- caret::trainControl(
      method = "repeatedcv",
      number = stab_params$tuning$folds,
      repeats = stab_params$tuning$repeats,
      summaryFunction = r_summary,
      allowParallel = FALSE
    )
    
    global_tune_grid <- NULL

        # Use the full training data set to generate a fixed lambda sequence. Then, lambda will be re-tuned in each iteration from the fixed sequence.
        glmnet_fit_global <- glmnet::glmnet(
          x = x_predictors, y = y_outcome, alpha = 1,
          nlambda = stab_params$tuning$tuneLength
        )
        global_tune_grid <- expand.grid(alpha = 1, lambda = glmnet_fit_global$lambda)

    n_cores <- config$main_pipeline_params$n_cores
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)

    set.seed(config$main_pipeline_params$seed)
    lgr::lgr$info(paste("Running", stab_params$iterations, "iterations..."))

    all_selected_probes <- foreach::foreach(
        i = 1:stab_params$iterations,
        .combine = 'c',
        .packages = c('glmnet', 'caret')
    ) %dorng% {
      sub_indices <- sample(1:nrow(x_predictors), size = floor(stab_params$fraction * nrow(x_predictors)))
      x_sub <- x_predictors[sub_indices, , drop = FALSE]
      y_sub <- y_outcome[sub_indices]

      final_model_iter <- NULL
      
        # Tune on the subsample to find the best lambda from the fixed sequence
        tune_grid_for_iter <- global_tune_grid
        if (is.null(tune_grid_for_iter)) {
            glmnet_fit_iter <- glmnet::glmnet(
                x = x_sub, y = y_sub, alpha = 1,
                nlambda = stab_params$tuning$tuneLength
            )
            tune_grid_for_iter <- expand.grid(alpha = 1, lambda = glmnet_fit_iter$lambda)
        }
      
        final_model_iter <- caret::train(
            x = x_sub, y = y_sub, method = "glmnet",
            trControl = train_control_tuning,
            tuneGrid = tune_grid_for_iter,
            metric = lambda_metric,
            verbose = FALSE
        )
      
      # Extract coefficients from the final model
      coefs <- coef(final_model_iter$finalModel, s = final_model_iter$bestTune$lambda)
      selected_probes <- rownames(coefs)[which(coefs[, 1] != 0)]
    }

    parallel::stopCluster(cl)
    foreach::registerDoSEQ()

  # Get selection probabilities
    if (length(all_selected_probes) > 0) {
      all_selected_probes <- all_selected_probes[all_selected_probes != "(Intercept)"]
      counts <- table(all_selected_probes)
      selection_prob <- setNames(rep(0, ncol(x_predictors)), colnames(x_predictors))
      selection_prob[names(counts)] <- counts / stab_params$iterations
    } else {
      selection_prob <- setNames(rep(0, ncol(x_predictors)), colnames(x_predictors))
    }

    stability_df <- data.frame(probe = names(selection_prob), raw_importance = as.numeric(selection_prob))
    save_ranking(stability_df, "ranking_lasso_stability_selection.csv", rank_probes_dir, "lasso_stability_selection")
}

  lgr::lgr$info("Probe ranking complete.")
  return(rank_probes_dir)
}