# apply_ensemble.R
# This function can be used to generate predictions from the cortisol response methylation profile
# score (MPS) described in 'Development and cross-tissue validation of a machine learning-based
# methylation profile score for the cortisol response to stress' (Balfour et al., 2026).

# new_data_path: Path to the new data, as an RDS file. Should be a list containing:
#                'pheno':  A data frame that includes the target variable (sample ID should be in
#                          the row names and phenotype variable names should be in the column
#                          names)
#                'methyl': A matrix of beta values or M values (sample ID should be in the column
#                          names and probe IDs should be in the row names) - must be from the
#                          Illumina Infinium 450K or EPIC BeadChip
# model_path: Path to the ensemble RDS file (cortisol_response_mps_2026.rds)
# impute: String or NULL specifying how to handle missing probes. Can be:
#         NULL: Missing probes are passed to the model as NAs
#         "blood": Imputes missing probes using pre-calculated medians from GSE77445
#         "saliva": Imputes missing probes using pre-calculated medians from GSE119078
#         File path: Path to an RDS list containing a 'methyl' beta matrix for imputation
# beta_to_m: If methylation values should be converted to M values. Set to TRUE if the object at
#            new_data_path contains beta values.
# remove_methyl_outliers: If TRUE, temporarily removes 'rs' probes and probes with missingness >=
#                         .05 and then removes samples where methylation PC1 or PC2 z score is > 3.
# samples_to_remove: Vector of characters specifying any additional sample IDs to remove.
# target_pheno: A vector of character strings specifying the target variables in 'pheno'. If
#               provided, applies winsorization (1.5 IQR threshold) and z-score conversion to all
#               specified target variables, saving them to the returned object with a '_proc'
#               suffix. Also runs a correlation test between these variables and model predictions.
# results_dir: Path to directory where the correlation results CSV file will be saved.

# This function returns 'pheno' from the object at new_data_path with the base model predictions,
# ensemble predictions, and - if target_pheno is provided - the processed target variables added.
# Also performs correlation tests between model predictions and the processed variables in
# target_pheno. The results are saved to a CSV file in results_dir. Please note that this function
# requires caret (version 6.0.94) and xgboost (version 1.7.1.1).

# Example:
# results <- apply_ensemble(
#    new_data_path = "./champagne_pre_formatted_beta.rds",
#    model_path = "./cortisol_response_mps_2026.rds",
#    impute = "saliva",
#    beta_to_m = TRUE,
#    remove_methyl_outliers = TRUE,
#    samples_to_remove = c("205851990031_R07C01", "205837550052_R03C01"),
#    target_pheno = "AUCi_scc",
#    results_dir = "./results"
# )

apply_ensemble <- function(new_data_path, model_path = "./cortisol_response_mps_2026.rds",
                           impute = NULL, beta_to_m = TRUE, remove_methyl_outliers = TRUE,
                           samples_to_remove = NULL, target_pheno = NULL, results_dir = NULL) {
  
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("The package 'caret' (version 6.0.94) is required but not installed.")
  } else if (packageVersion("caret") != "6.0.94") {
    stop(sprintf("Error: 'caret' version 6.0.94 is required, but you have version %s.", packageVersion("caret")))
  }

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("The package 'xgboost' (version 1.7.1.1) is required but not installed.")
  } else if (packageVersion("xgboost") != "1.7.1.1") {
    stop(sprintf("Error: 'xgboost' version 1.7.1.1 is required, but you have version %s.", packageVersion("xgboost")))
  }  
  
  suppressPackageStartupMessages({
    library(caret)
    library(xgboost)
  })

  message("Validating inputs...")  
  
  if (!is.null(target_pheno) && is.null(results_dir)) {
    stop("Error: 'target_pheno' was provided but 'results_dir' is NULL. You must specify a directory to save the correlations.")
  }
  
  new_data <- readRDS(new_data_path)
  
  if (!is.list(new_data) || !all(c("pheno", "methyl") %in% names(new_data))) {
    stop("Error: 'new_data' must be a list containing 'pheno' and 'methyl'.")
  }
  
  pheno <- new_data$pheno
  methyl <- new_data$methyl
  
  if (!identical(colnames(methyl), rownames(pheno))) {
    stop("Error: Columns of 'methyl' do not match rows of 'pheno'.")
  }
  
  if (!is.null(target_pheno)) {
    missing_phenos <- setdiff(target_pheno, colnames(pheno))
    if (length(missing_phenos) > 0) {
      stop(sprintf("Error: Target phenotype/s not found in pheno data frame: %s", paste(missing_phenos, collapse = ", ")))
    }
  }
  
  if (beta_to_m) {
    min_val <- min(methyl, na.rm = TRUE)
    max_val <- max(methyl, na.rm = TRUE)
    if (min_val < 0 || max_val > 1) {
      stop(sprintf("Error: Target methylation data not within the expected range for beta values (range 0 to 1). Min: %.4f, Max: %.4f", min_val, max_val))
    }
  }
  
  message("Loading ensemble model...")
  if (!file.exists(model_path)) stop("Error: Model not found at ", model_path)
  model_object <- suppressWarnings(readRDS(model_path))
  
  # Strict model validation check
  if (is.null(model_object$mps_signature) || model_object$mps_signature != "cortisol_response_mps_2026") {
    stop("Error: The provided object is not the valid cortisol_response_mps_2026.rds model.")
  }
  
  base_models <- model_object$base_models
  meta_model <- model_object$meta_model
  
  required_probes <- unique(unlist(lapply(base_models, function(m) m$model$coefnames)))
  
  if (remove_methyl_outliers) {
    message("Running PCA-based methylation outlier detection...")
    methyl_for_pca <- methyl
    
    prefixes <- "rs"
    probe_miss_thresh <- 0.05
    threshold <- 3
    
    if (length(prefixes) > 0) {
      pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
      probes_to_flag <- rownames(methyl_for_pca)[grepl(pattern, rownames(methyl_for_pca))]
      probes_to_remove_pca <- setdiff(probes_to_flag, required_probes)
      if (length(probes_to_remove_pca) > 0) methyl_for_pca <- methyl_for_pca[!rownames(methyl_for_pca) %in% probes_to_remove_pca, , drop = FALSE]
    }
    
    if (!is.null(probe_miss_thresh)) {
      row_miss <- rowMeans(is.na(methyl_for_pca))
      methyl_for_pca <- methyl_for_pca[row_miss < probe_miss_thresh, , drop = FALSE]
    }
    
    if (any(is.na(methyl_for_pca))) {
      methyl_for_pca <- apply(methyl_for_pca, 2, function(col) {
        col[is.na(col)] <- mean(col, na.rm = TRUE)
        col
      })
    }
    
    row_means <- rowMeans(methyl_for_pca, na.rm = TRUE)
    row_sq_dev <- (methyl_for_pca - row_means)^2 
    n_obs <- rowSums(!is.na(methyl_for_pca))
    n_obs[n_obs <= 1] <- NA
    probe_variances <- rowSums(row_sq_dev, na.rm = TRUE) / (n_obs - 1)
    
    probes_to_keep_pca <- !is.na(probe_variances) & probe_variances > 1e-8 
    methyl_for_pca <- methyl_for_pca[probes_to_keep_pca, , drop = FALSE]
    
    pca_result <- prcomp(t(methyl_for_pca), scale. = TRUE, center = TRUE)
    
    pc1_zscore <- as.numeric(scale(pca_result$x[, 1]))
    pc2_zscore <- as.numeric(scale(pca_result$x[, 2]))
    outlier_idx <- which(abs(pc1_zscore) > threshold | abs(pc2_zscore) > threshold)
    
    if (length(outlier_idx) > 0) {
      outlier_samples <- colnames(methyl)[outlier_idx]
      message(sprintf("Removed %d methylation outlier sample(s).", length(outlier_samples)))
      samps_to_keep <- setdiff(colnames(methyl), outlier_samples)
      methyl <- methyl[, samps_to_keep, drop = FALSE]
      pheno <- pheno[samps_to_keep, , drop = FALSE]
    } else {
      message("No methylation outlier samples found.")
    }
  }
  
  if (!is.null(samples_to_remove) && length(samples_to_remove) > 0) {
    current_samples <- colnames(methyl)
    samples_to_manually_remove <- intersect(samples_to_remove, current_samples)
    if (length(samples_to_manually_remove) > 0) {
      message(sprintf("Removing %d manually specified sample/s.", length(samples_to_manually_remove)))
      samps_to_keep <- setdiff(current_samples, samples_to_manually_remove)
      methyl <- methyl[, samps_to_keep, drop = FALSE]
      pheno <- pheno[samps_to_keep, , drop = FALSE]
    }
  }
  
  present_probes <- rownames(methyl)
  missing_probes <- setdiff(required_probes, present_probes)
  n_samples <- ncol(methyl)
  
  if (length(missing_probes) > 0) {
    message(sprintf("Found %d missing required probes.", length(missing_probes)))
    
    if (!is.null(impute)) {
      if (impute %in% c("blood", "saliva")) {
        message(sprintf("Using built-in %s imputation medians...", impute))
        ref_values <- model_object$imputation_medians[[impute]]
        
        if (is.null(ref_values)) {
          stop(sprintf("Error: %s imputation medians not found in the model object.", impute))
        }
        
        global_fallback_value <- median(ref_values, na.rm = TRUE)
        
      } else {
        message("Using provided imputation reference dataset...")
        if (!file.exists(impute)) stop("Error: Imputation reference file not found at ", impute)
        
        ref_data <- readRDS(impute)
        
        if (beta_to_m) {
          ref_min <- min(ref_data$methyl, na.rm = TRUE)
          ref_max <- max(ref_data$methyl, na.rm = TRUE)
          if (ref_min < 0 || ref_max > 1) {
            stop(sprintf("Error: beta_to_m is TRUE, but imputation reference data appears to be M values. Both target and reference data must be beta values. Reference min: %.4f, max: %.4f", ref_min, ref_max))
          }
        }
        
        available_in_ref <- intersect(missing_probes, rownames(ref_data$methyl))
        if (length(available_in_ref) > 0) {
          ref_values <- apply(ref_data$methyl[available_in_ref, , drop = FALSE], 1, median, na.rm = TRUE)
        } else {
          ref_values <- numeric(0)
        }
        
        global_pool_probes <- intersect(required_probes, rownames(ref_data$methyl))
        if (length(global_pool_probes) == 0) stop("Error: Imputation reference contains no required probes.")
        global_fallback_value <- median(as.numeric(ref_data$methyl[global_pool_probes, , drop = FALSE]), na.rm = TRUE)
      }
      
      available_in_ref_values <- intersect(missing_probes, names(ref_values))
      still_missing <- setdiff(missing_probes, available_in_ref_values)
      
      imputed_chunks <- list()
      
      if (length(available_in_ref_values) > 0) {
        chunk_matrix <- matrix(rep(ref_values[available_in_ref_values], times = n_samples), nrow = length(available_in_ref_values), ncol = n_samples, byrow = FALSE)
        rownames(chunk_matrix) <- available_in_ref_values
        imputed_chunks[[1]] <- chunk_matrix
      }
      
      if (length(still_missing) > 0) {
        chunk_matrix_fallback <- matrix(global_fallback_value, nrow = length(still_missing), ncol = n_samples)
        rownames(chunk_matrix_fallback) <- still_missing
        imputed_chunks[[2]] <- chunk_matrix_fallback
      }
      
      full_imputed_matrix <- do.call(rbind, imputed_chunks)
      colnames(full_imputed_matrix) <- colnames(methyl)
      methyl <- rbind(methyl, full_imputed_matrix)
      
    } else {
      message("No imputation reference provided.")
      na_matrix <- matrix(NA, nrow = length(missing_probes), ncol = n_samples)
      rownames(na_matrix) <- missing_probes
      colnames(na_matrix) <- colnames(methyl)
      methyl <- rbind(methyl, na_matrix)
    }
  }
  
  methyl_req <- methyl[required_probes, , drop = FALSE]
  
  if (beta_to_m) {
    message("Transforming beta values to M values for required probes...")
    offset <- 0.001
    
    methyl_req[methyl_req == 0] <- offset
    methyl_req[methyl_req == 1] <- 1 - offset
    methyl_req <- log2(methyl_req / (1 - methyl_req))
  }
  
  if (!is.null(target_pheno)) {
    message("Processing target phenotype(s)...")
    for (tp in target_pheno) {
      pheno_target <- pheno[[tp]]
      
      # Hardcoded processing parameters
      n_iqr <- 1.5
      do_stand <- TRUE
      
      q1 <- unname(quantile(pheno_target, 0.25, na.rm = TRUE))
      q3 <- unname(quantile(pheno_target, 0.75, na.rm = TRUE))
      iqr <- q3 - q1
      lower_bound <- q1 - n_iqr * iqr
      upper_bound <- q3 + n_iqr * iqr
      
      pheno_target[which(pheno_target < lower_bound)] <- lower_bound
      pheno_target[which(pheno_target > upper_bound)] <- upper_bound
      
      if (do_stand) {
        pheno_df_temp <- data.frame(pheno = pheno_target)
        stand_model <- caret::preProcess(pheno_df_temp, method = c("center", "scale"))
        pheno_target <- predict(stand_model, newdata = pheno_df_temp)$pheno
      }
      
      pheno[[paste0(tp, "_proc")]] <- pheno_target
    }
  }
  
  message("Generating predictions from base models...")
  
  base_predictions_list <- lapply(seq_along(base_models), function(i) {
    bm <- base_models[[i]]$model
    bm_name <- names(base_models)[i]
    if (is.null(bm_name) || bm_name == "") bm_name <- paste0("Model_", i)
    
    bm_probes <- bm$coefnames    
    pred_matrix <- t(methyl_req[bm_probes, , drop = FALSE])
    
    fm <- bm$finalModel

    # Reconstruct booster
    native_booster <- suppressWarnings(xgboost::xgb.load.raw(fm))
    preds <- predict(native_booster, newdata = pred_matrix)
    
    return(preds)
  })
  
  base_predictions_df <- as.data.frame(base_predictions_list)
  
  model_names <- names(base_models)
  if (is.null(model_names)) model_names <- paste0("BaseModel_", seq_along(base_models))
  
  colnames(base_predictions_df) <- model_names 
  rownames(base_predictions_df) <- colnames(methyl_req)
  
  message("Generating ensemble predictions...")
  final_predictions <- predict(meta_model, newdata = base_predictions_df, na.action = na.pass)
  
  colnames(base_predictions_df) <- paste0(model_names, "_pred")
  
  pheno <- cbind(pheno, base_predictions_df)
  pheno$ensemble_pred <- final_predictions
  
  if (!is.null(target_pheno)) {
    message("Calculating correlations...")
    cor_list <- list()
    pred_cols <- c(paste0(model_names, "_pred"), "ensemble_pred")
    
    for (tp in target_pheno) {
      tp_proc <- paste0(tp, "_proc")
      for (pcol in pred_cols) {
        # Select complete cases for testing
        valid <- complete.cases(pheno[[tp_proc]], pheno[[pcol]])
        
        if (sum(valid) > 2) {
          ctest <- cor.test(pheno[[tp_proc]][valid], pheno[[pcol]][valid])
          cor_list[[length(cor_list) + 1]] <- data.frame(
            Model = pcol,
            Target_Variable = tp,
            Estimate = ctest$estimate,
            P_Value = ctest$p.value,
            stringsAsFactors = FALSE
          )
        } else {
          cor_list[[length(cor_list) + 1]] <- data.frame(
            Model = pcol,
            Target_Variable = tp,
            Estimate = NA,
            P_Value = NA,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    
    cor_df <- do.call(rbind, cor_list)
    cor_df$FDR_P_Value <- p.adjust(cor_df$P_Value, method = "fdr")    
    cor_df <- cor_df[order(cor_df$FDR_P_Value), ]
    
    if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)
    out_path <- file.path(results_dir, "cortisol_response_mps_correlations.csv")
    write.csv(cor_df, out_path, row.names = FALSE)
    message(sprintf("Saved correlation results to %s", out_path))
  }
  
  message("Done.")
  return(invisible(pheno))
}