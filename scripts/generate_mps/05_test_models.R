# 05_test_models.R
test_models <- function(config, data, train_models_dir, out_dir) {
  
  lgr::lgr$info("Starting model testing...")
  
  test_models_dir <- file.path(out_dir, "05_test_models")
  if (!dir.exists(test_models_dir)) dir.create(test_models_dir, recursive = TRUE)
  
  scatter_plots_dir <- file.path(test_models_dir, "scatter_plots")
  if (!dir.exists(scatter_plots_dir)) dir.create(scatter_plots_dir, recursive = TRUE)

  processed_test_data <- data  
  model_files <- list.files(train_models_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  
  # Loop through models, predict, and evaluate
  results_list <- purrr::map(model_files, function(model_file) {
    lgr::lgr$info(paste("Testing model:", basename(model_file)))
    
    # Load model object (including metadata from train_models())
    model_obj <- readRDS(model_file)
    model <- model_obj$model
    meta <- model_obj$metadata
    
    full_method_name <- meta$ranking_method
    
    probes_before_filtering_match <- regexpr("_[0-9]+", full_method_name)
    probes_before_filtering <- if (probes_before_filtering_match != -1) {
      as.integer(gsub("_", "", regmatches(full_method_name, probes_before_filtering_match)))
    } else {
      meta$n_probes 
    }
    
    base_ranking_method <- gsub("_[0-9]+.*$", "", full_method_name)
    
    probes_in_model <- model$coefnames
    
# Get samples that are present in both pheno and methyl
common_samples <- intersect(rownames(processed_test_data$pheno), 
                            colnames(processed_test_data$methyl))

# Filter both data sets to the common samples
pheno_aligned <- processed_test_data$pheno[common_samples, , drop = FALSE]
methyl_subset <- t(processed_test_data$methyl[probes_in_model, common_samples, drop = FALSE])

# Check if the row names still match, to ensure data alignment
if (!all(rownames(pheno_aligned) == rownames(methyl_subset))) {
  stop("Sample names do not align.")
}

test_df <- data.frame(
  pheno = pheno_aligned[[config$process_data_params$pheno$pheno_name]],
  methyl_subset
)

    predictions <- predict(model, newdata = test_df)
eval_data <- data.frame(obs = test_df$pheno, pred = predictions)

# Calculate metrics
all_metrics <- r_summary(eval_data)

results_df <- data.frame(
  path = model_file,
  outer_fold = meta$outer_fold,
  ranking_method = base_ranking_method,
  training_method = meta$training_method,
  n_probes = probes_before_filtering,
  Rsquared = all_metrics["Rsquared"],
  Rsquared_ss = all_metrics["Rsquared_ss"],
  r = all_metrics["r"],
  rho = all_metrics["rho"],
  RMSE = all_metrics["RMSE"],
  MAE = all_metrics["MAE"]
)    

    plot_df <- data.frame(
      observed = test_df$pheno,
      predicted_orig = predictions
    )
    
    model_basename <- tools::file_path_sans_ext(basename(model_file))
    
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(y = observed)) +
      ggplot2::geom_point(ggplot2::aes(x = predicted_orig), color = "black", alpha = 0.7) +
      ggplot2::geom_smooth(ggplot2::aes(x = predicted_orig), method = "lm", se = FALSE, color = "black", formula = y ~ x) +
      ggplot2::labs(
        title = paste("Observed vs. predicted values for", model_basename),
        subtitle = paste("Ranking:", base_ranking_method, "| Training:", meta$training_method, "| Probes:", probes_before_filtering),
        x = "Predicted",
        y = "Observed"
      ) +
ggplot2::annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -1.5,
                  label = paste("Original r =", round(all_metrics["r"], 3)), color = "black") +
      ggplot2::theme_minimal() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5))
      
    plot_filename <- file.path(scatter_plots_dir, paste0(model_basename, "_scatter.pdf"))
    ggplot2::ggsave(plot_filename, p, width = 8, height = 6)
    
    return(results_df)
  })
  
  final_results_df <- dplyr::bind_rows(results_list)
  
  utils::write.csv(final_results_df, file.path(test_models_dir, "test_results.csv"), row.names = FALSE)
  
  lgr::lgr$info("Model testing complete.")
  return(test_models_dir)
}