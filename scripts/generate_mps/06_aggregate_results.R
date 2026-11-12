# 06_aggregate_results.R
aggregate_results <- function(config, out_dir) {
  
  lgr::lgr$info("Starting results aggregation...")
  
  result_files <- list.files(
    path = out_dir, pattern = "test_results.csv",
    recursive = TRUE, full.names = TRUE
  )
  
  if (length(result_files) == 0) {
    lgr::lgr$warn("No test result files found to aggregate.")
    return(NULL)
  }
  
  all_results <- purrr::map_dfr(result_files, utils::read.csv)
  
  performance_dir <- file.path(out_dir, "performance")
  if (!dir.exists(performance_dir)) dir.create(performance_dir)
  
  # Save sorted results
  if (is.null(config$aggregate_results_params$sorting_metric)) {
    stop("sorting_metric must be specified in aggregate_results_params")
  }
  sorting_metric <- config$aggregate_results_params$sorting_metric
  lgr::lgr$info(paste("Sorting aggregated results by:", sorting_metric))
  
  maximise <- config$metric_maximise_lookup[[sorting_metric]]
  if (is.null(maximise)) stop(paste("Metric", sorting_metric, "not found in config."))
  
  if (maximise) {
    sorted_results <- all_results %>% dplyr::arrange(dplyr::desc(.data[[sorting_metric]]))
  } else {
    sorted_results <- all_results %>% dplyr::arrange(.data[[sorting_metric]])
  }
  
  main_output_path <- file.path(performance_dir, "aggregated_test_results.csv")
  utils::write.csv(sorted_results, main_output_path, row.names = FALSE)
  lgr::lgr$info(paste("Aggregated test results saved to:", main_output_path))
  
  lgr::lgr$info("Results aggregation complete.")
  return(main_output_path)
}