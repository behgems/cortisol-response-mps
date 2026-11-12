# 01_process_data.R
process_data <- function(config, data, out_dir, process_methyl = TRUE, process_pheno = TRUE, required_probes = NULL, sub_dir_name = "01_process_data", samples_to_remove = NULL) {
  
  lgr::lgr$info("Starting data processing...")
  
  process_data_dir <- file.path(out_dir, sub_dir_name)
  if (!dir.exists(process_data_dir)) dir.create(process_data_dir, recursive = TRUE)
  
  pheno_params <- config$process_data_params$pheno
  methyl_params <- config$process_data_params$methyl
  
  pheno_df <- data$pheno
  methyl_mat <- data$methyl
  
  pheno_name <- pheno_params$pheno_name
  
  if (process_methyl) {
    lgr::lgr$info("Starting methylation processing...")
    
    probes_to_flag_for_removal <- c()

    # Prefix probe filtering
    if (!is.null(methyl_params$filter_probe_prefix) && length(methyl_params$filter_probe_prefix) > 0) {
      prefixes <- methyl_params$filter_probe_prefix
      pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
      probes_from_prefix_filter <- rownames(methyl_mat)[grepl(pattern, rownames(methyl_mat))]
      probes_to_flag_for_removal <- c(probes_to_flag_for_removal, probes_from_prefix_filter)
      lgr::lgr$info(paste0("Flagged ", length(probes_from_prefix_filter), " probes for removal by prefix(es): '", 
                           paste(prefixes, collapse = ", "), "'."))
    }

    probes_to_flag_for_removal <- unique(probes_to_flag_for_removal)
    
    if (length(probes_to_flag_for_removal) > 0) {
      probes_to_remove <- setdiff(probes_to_flag_for_removal, required_probes) # Protect probes required by model
      
      if (length(probes_to_remove) > 0) {
        initial_probes <- nrow(methyl_mat)
        methyl_mat <- methyl_mat[!rownames(methyl_mat) %in% probes_to_remove, ]
        lgr::lgr$info(paste0("Removed a total of ", initial_probes - nrow(methyl_mat), " probes."))
      } else {
         lgr::lgr$info("No probes were removed.")
      }
    } else {
        lgr::lgr$info("No probe filtering requested.")
    }
    
    # Filter samples by missingness
    probe_miss_thresh <- methyl_params$probe_miss_thresh
    if (!is.null(probe_miss_thresh)) {
      lgr::lgr$info(paste0("Filtering probes with missingness > ", probe_miss_thresh, "."))
      initial_probes <- nrow(methyl_mat)
      row_miss <- rowMeans(is.na(methyl_mat))
      
      probes_to_keep <- names(row_miss[row_miss < probe_miss_thresh])
      
      methyl_mat <- methyl_mat[probes_to_keep, , drop = FALSE]
      lgr::lgr$info(paste0("Removed ", initial_probes - nrow(methyl_mat), " probe(s). ", nrow(methyl_mat), " probes retained."))
    } else {
      lgr::lgr$info("Probe missingness threshold is NULL, so skipping probe filtering.")
    }
    
    # PCA-based sample outlier handling
    if (!is.null(methyl_params$sample_outlier_handling) && isTRUE(methyl_params$sample_outlier_handling$use)) {
      lgr::lgr$info("Starting PCA-based methylation outlier handling...")
      threshold <- methyl_params$sample_outlier_handling$threshold
      lgr::lgr$info(paste("Using a z-score threshold of", threshold, "on PC1 or PC2..."))
      
      if (any(is.na(methyl_mat))) {
        lgr::lgr$info("Performing temporary mean imputation for PCA...") # Cannot have NAs when doing PCA
        methyl_mat_for_pca <- apply(methyl_mat, 2, function(col) {
          col[is.na(col)] <- mean(col, na.rm = TRUE)
          col
        })
      } else {
        methyl_mat_for_pca <- methyl_mat
      }
      
      lgr::lgr$info("Checking for zero-variance probes before running PCA...")
      
      row_means <- rowMeans(methyl_mat_for_pca, na.rm = TRUE)
      row_sq_dev <- (methyl_mat_for_pca - row_means)^2 
      n_obs <- rowSums(!is.na(methyl_mat_for_pca))
      n_obs[n_obs <= 1] <- NA
      probe_variances <- rowSums(row_sq_dev, na.rm = TRUE) / (n_obs - 1)
      
      probes_to_keep_pca <- !is.na(probe_variances) & probe_variances > 1e-8 
      
      n_removed <- sum(!probes_to_keep_pca)
      if (n_removed > 0) {
        lgr::lgr$info(paste("Found and temporarily removed", n_removed, "zero-variance probe(s) for PCA."))
        methyl_mat_for_pca <- methyl_mat_for_pca[probes_to_keep_pca, ]
      } else {
        lgr::lgr$info("No zero-variance probes found.")
      }
      
      lgr::lgr$info("Running PCA on methylation data...")
      pca_result <- prcomp(t(methyl_mat_for_pca), scale. = TRUE, center = TRUE)
      pca_data <- as.data.frame(pca_result$x[, 1:2])
      
      pca_data$pc1_zscore <- as.numeric(scale(pca_data$PC1))
      pca_data$pc2_zscore <- as.numeric(scale(pca_data$PC2))
      outlier_idx <- which(abs(pca_data$pc1_zscore) > threshold | abs(pca_data$pc2_zscore) > threshold)
      
      pca_data$sample_id <- colnames(methyl_mat)
      pca_data$is_outlier <- ifelse(1:nrow(pca_data) %in% outlier_idx, "yes", "no")
      
      variance_explained <- summary(pca_result)$importance[2, 1:2] * 100
      
      pca_plot <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = is_outlier, label = sample_id)) +
        ggplot2::geom_point(alpha = 0.8) +
        ggrepel::geom_text_repel(data = subset(pca_data, is_outlier == "yes"), max.overlaps = Inf, box.padding = 0.5) +
        ggplot2::scale_color_manual(values = c("no" = "black", "yes" = "red"), name = "Outlier") +
        ggplot2::labs(
          title = "Methylation PCA for outlier detection",
          subtitle = paste("Outliers identified using a", threshold, "z-score threshold on PC1 or PC2"),
          x = sprintf("PC 1 (%.2f%%)", variance_explained[1]),
          y = sprintf("PC 2 (%.2f%%)", variance_explained[2])
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")
      
      pdf_path <- file.path(process_data_dir, "methyl_pca_outlier_plot.pdf")
      ggplot2::ggsave(filename = pdf_path, plot = pca_plot, width = 8, height = 7)
      lgr::lgr$info(paste("Saved PCA outlier plot to:", pdf_path))
      
      if (length(outlier_idx) > 0) {
        outlier_samples <- colnames(methyl_mat)[outlier_idx]
        lgr::lgr$info(paste0("Identified ", length(outlier_samples), " outlier sample(s): ", paste(outlier_samples, collapse = ", ")))
        
        outliers_df <- pca_data[outlier_idx, c("sample_id", "pc1_zscore", "pc2_zscore")]
        outliers_df$distance_from_mean <- sqrt(outliers_df$pc1_zscore^2 + outliers_df$pc2_zscore^2)
        
        csv_path <- file.path(process_data_dir, "methyl_outliers_report.csv")
        write.csv(outliers_df, file = csv_path, row.names = FALSE)
        lgr::lgr$info(paste("Saved methylation outlier report to:", csv_path))
        
        initial_samps <- ncol(methyl_mat)
        samps_to_keep <- setdiff(colnames(methyl_mat), outlier_samples)
        methyl_mat <- methyl_mat[, samps_to_keep, drop = FALSE]
        if (!is.null(pheno_df)) {
          pheno_df <- pheno_df[samps_to_keep, , drop = FALSE]
        }
        lgr::lgr$info(paste0("Removed ", initial_samps - length(samps_to_keep), " outlier sample(s). ", length(samps_to_keep), " samples retained."))
        
      } else {
        lgr::lgr$info("No methylation outliers found.")
      }
    }
   
    # Manual sample removal (to handle situation where sample is an outlier in post but not pre or vice-versa)
    if (!is.null(samples_to_remove) && length(samples_to_remove) > 0) {
      lgr::lgr$info("Starting manual sample removal...")
      
      # Check which manually specified samples are currently in the dataset
      current_samples <- colnames(methyl_mat)
      samples_to_manually_remove <- intersect(samples_to_remove, current_samples)
      
      # Log samples that were requested but not found
      samples_not_found <- setdiff(samples_to_remove, current_samples)
      if (length(samples_not_found) > 0) {
        lgr::lgr$warn(paste0("The following manually specified samples were not found: ", 
                             paste(samples_not_found, collapse = ", ")))
      }
      
      if (length(samples_to_manually_remove) > 0) {
        lgr::lgr$info(paste0("Removing ", length(samples_to_manually_remove), " specified sample(s): ", 
                             paste(samples_to_manually_remove, collapse = ", ")))
        
        initial_samps <- ncol(methyl_mat)
        samps_to_keep <- setdiff(current_samples, samples_to_manually_remove)
        
        # Filter methyl matrix
        methyl_mat <- methyl_mat[, samps_to_keep, drop = FALSE]
        
        # Filter pheno df
        if (!is.null(pheno_df)) {
          pheno_df <- pheno_df[samps_to_keep, , drop = FALSE]
        }
        
        lgr::lgr$info(paste0("Removed ", initial_samps - length(samps_to_keep), " manually specified sample(s). ", length(samps_to_keep), " samples retained."))
        
      } else {
        lgr::lgr$info("No samples to remove from the manual list (all either not found or already removed).")
      }
    } else {
      if (process_methyl) {
         lgr::lgr$info("No manual sample removal requested (i.e., 'samples_to_remove' is NULL or empty).")
      }
    }

    # Betas to M values
    if (methyl_params$beta_to_m$use) {
      min_val <- min(methyl_mat, na.rm = TRUE)
      max_val <- max(methyl_mat, na.rm = TRUE)
      if (min_val < 0 || max_val > 1) {
        lgr::lgr$warn(paste0("Beta to M conversion is enabled, but values were found outside the expected range for betas (0 to 1). Min: ", round(min_val, 4), ", Max: ", round(max_val, 4), "."))
      }
      
      lgr::lgr$info("Converting betas to M values...")
        offset <- methyl_params$beta_to_m$offset
        methyl_mat[methyl_mat == 0] <- offset
        methyl_mat[methyl_mat == 1] <- 1 - offset
        methyl_mat <- log2(methyl_mat / (1 - methyl_mat))
      
    }
  }
  
  # Phenotype processing
  if (process_pheno) {
    lgr::lgr$info("Starting phenotype processing...")
    
    if (is.null(pheno_df) || !pheno_name %in% colnames(pheno_df)) {
      stop(paste("Phenotype '", pheno_name, "' not found in pheno data frame."))
    }
    pheno_target <- pheno_df[[pheno_name]]
    
    if (pheno_params$outlier_handling$use) {
      method <- pheno_params$outlier_handling$method
      lgr::lgr$info(paste0("Handling phenotype outliers using winsorisation with IQR of ", pheno_params$outlier_handling$n_iqr))

      q1 <- unname(quantile(pheno_target, 0.25, na.rm = TRUE))
      q3 <- unname(quantile(pheno_target, 0.75, na.rm = TRUE))
      iqr <- q3 - q1
      lower_bound <- q1 - pheno_params$outlier_handling$n_iqr * iqr
      upper_bound <- q3 + pheno_params$outlier_handling$n_iqr * iqr
      outliers_idx <- which(pheno_target < lower_bound | pheno_target > upper_bound)
      
      plot_df <- data.frame(sample_id = rownames(pheno_df), pheno_value = pheno_target)
      plot_df$label <- ""
      if (length(outliers_idx) > 0) {
        plot_df$label[outliers_idx] <- plot_df$sample_id[outliers_idx]
      }
      
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = "", y = pheno_value)) +
          ggplot2::geom_boxplot(outlier.shape = NA) +
          ggplot2::geom_jitter(color = "blue", alpha = 0.6, width = 0.1) +
          ggrepel::geom_text_repel(ggplot2::aes(label = label), max.overlaps = Inf) +
          ggplot2::geom_hline(yintercept = upper_bound, linetype = "dashed", color = "red") +
          ggplot2::geom_hline(yintercept = lower_bound, linetype = "dashed", color = "red") +
          ggplot2::labs(
            title = "Distribution before outlier handling",
            subtitle = paste("Outliers identified using", pheno_params$outlier_handling$n_iqr, "IQR rule"),
            y = pheno_name,
            x = ""
          ) +
          ggplot2::theme_minimal()
      
      pdf_path <- file.path(process_data_dir, "pheno_dist_before_outlier_handling.pdf")
      ggplot2::ggsave(filename = pdf_path, plot = p, width = 8, height = 6)
      lgr::lgr$info(paste("Saved phenotype distribution boxplot to:", pdf_path))
      
      if (length(outliers_idx) > 0) {
        outlier_data <- lapply(outliers_idx, function(idx) {
          value <- pheno_target[idx]
          type <- ifelse(value > upper_bound, "upper", "lower")
          distance <- ifelse(type == "upper", (value - upper_bound) / iqr, (lower_bound - value) / iqr)
          data.frame(
            sample_id = rownames(pheno_df)[idx],
            raw_phenotype_value = value,
            outlier_type = type,
            distance_in_iqr = distance
          )
        })
        outliers_df <- do.call(rbind, outlier_data)
        csv_path <- file.path(process_data_dir, "pheno_outliers_report.csv")
        write.csv(outliers_df, file = csv_path, row.names = FALSE)
        lgr::lgr$info(paste("Saved outlier report to:", csv_path))
        
          outliers_lower_idx <- which(pheno_target < lower_bound)
          outliers_upper_idx <- which(pheno_target > upper_bound)
          lgr::lgr$info(paste0("Winsorising ", length(outliers_idx), " outlier(s)..."))
          pheno_target[outliers_lower_idx] <- lower_bound
          pheno_target[outliers_upper_idx] <- upper_bound
      } else {
        lgr::lgr$info("No outliers detected.")
      }
    }

   if (pheno_params$standardisation$use) {
      pheno_df_temp <- data.frame(pheno = pheno_target)
      stand_model <- caret::preProcess(pheno_df_temp, method = c("center", "scale"))
      pheno_target <- predict(stand_model, newdata = pheno_df_temp)$pheno
    }
    
    pheno_df[[pheno_name]] <- pheno_target
  }
  
  # Finalise and save
  lgr::lgr$info("Finalising processed data...")
  if (!is.null(pheno_df)) {
    stopifnot("Sample alignment mismatch" = all(rownames(pheno_df) == colnames(methyl_mat)))
  }
  
  processed_data <- list(pheno = pheno_df, methyl = methyl_mat)
  output_file_path <- file.path(process_data_dir, "processed_data.rds")
  saveRDS(processed_data, file = output_file_path)
  
  lgr::lgr$info(paste("Processed data saved to:", output_file_path))
  
  return(
    list(
      dir = process_data_dir
    )
  )
}