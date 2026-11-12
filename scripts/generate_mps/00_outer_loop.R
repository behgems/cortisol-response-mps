# 00_outer_loop.R
outer_loop <- function(config, data, out_dir) {
  
  if ("file" %in% names(lgr::lgr$appenders)) lgr::lgr$remove_appender("file")
  lgr::lgr$set_threshold("info")
  
  if (is.list(config$outer_loop_params$folds) && is.numeric(config$outer_loop_params$folds[[1]])) {
    original_sample_ids <- rownames(data$pheno)
    config$outer_loop_params$folds <- lapply(config$outer_loop_params$folds, function(fold_indices) {
      original_sample_ids[fold_indices]
    })
  }
  
  lgr::lgr$info(paste("Starting run at directory:", out_dir))
  dir.create(out_dir, recursive = TRUE)
  
  log_file <- file.path(out_dir, "log.txt")
  lgr::lgr$add_appender(lgr::AppenderFile$new(log_file), name = "file")
  lgr::lgr$info("===== 00 STARTING OUTER LOOP =====")
  
  set.seed(config$main_pipeline_params$seed)
  lgr::lgr$info(paste("Using seed:", config$main_pipeline_params$seed))
  lgr::lgr$info(paste("Using n_cores:", config$main_pipeline_params$n_cores))
  
  lgr::lgr$info("Processing methyl globally...")
  process_data_result <- process_data(config, data, out_dir, 
                                      process_methyl = TRUE, process_pheno = FALSE,
                                      sub_dir_name = "01_process_data/methyl")
  processed_data <- readRDS(file.path(process_data_result$dir, "processed_data.rds"))

  data_for_folds <- processed_data

  lgr::lgr$info("Using folds from config...")
  if (is.list(config$outer_loop_params$folds)) {
    lgr::lgr$info("Using list of folds...")
    folds <- config$outer_loop_params$folds 
  } else {
    stop("Invalid 'folds' parameter in config. Must be a list of sample IDs.")
  }

  lgr::lgr$info("Filtering folds to retain only samples present in the processed data...")
  initial_fold_samples <- unique(unlist(folds))
  samples_in_processed_data <- rownames(data_for_folds$pheno)
  
  folds <- lapply(folds, function(fold_samples) {
    intersect(fold_samples, samples_in_processed_data)
  })
  
  final_fold_samples <- unique(unlist(folds))
  lgr::lgr$info(paste("Removed", length(initial_fold_samples) - length(final_fold_samples), "sample(s) from fold definitions that were not in the processed data..."))
  
  lgr::lgr$info(paste("Running pipeline for all", length(folds), "folds..."))
  
  for (i in seq_along(folds)) {
    fold_name <- paste0("fold_", i)
    out_dir_fold <- file.path(out_dir, fold_name)
    if (!dir.exists(out_dir_fold)) dir.create(out_dir_fold)

    lgr::lgr$info("===== 01 PROCESSING DATA =====")
    lgr::lgr$info(paste("Processing pheno for fold", i))
      
    test_sample_ids <- folds[[i]]
    train_sample_ids <- setdiff(rownames(processed_data$pheno), test_sample_ids)
      
    lgr::lgr$info("Processing training set pheno...")
    train_data_for_pheno_proc <- list(pheno = processed_data$pheno[train_sample_ids, , drop = FALSE], methyl = processed_data$methyl[, train_sample_ids, drop = FALSE])
      
    train_pheno_result <- process_data(config, train_data_for_pheno_proc, out_dir_fold, 
                                       process_methyl = FALSE, process_pheno = TRUE,
                                       sub_dir_name = "01_process_data/train")
      
    lgr::lgr$info("Processing testing set pheno...")
    test_data_for_pheno_proc <- list(pheno = processed_data$pheno[test_sample_ids, , drop = FALSE], methyl = processed_data$methyl[, test_sample_ids, drop = FALSE])
      
    test_pheno_result <- process_data(config, test_data_for_pheno_proc, out_dir_fold, 
                                      process_methyl = FALSE, process_pheno = TRUE,
                                      sub_dir_name = "01_process_data/test")

    training_data_pipeline <- readRDS(file.path(train_pheno_result$dir, "processed_data.rds"))
    testing_data_pipeline <- readRDS(file.path(test_pheno_result$dir, "processed_data.rds"))
    
    fold_specific_config <- config
    fold_specific_config$main_pipeline_params$seed <- config$main_pipeline_params$seed + i
    
    lgr::lgr$info(paste("Using seed", fold_specific_config$main_pipeline_params$seed, "for fold", i))

    main_pipeline(config = fold_specific_config,
                  training_data = training_data_pipeline,
                  testing_data = testing_data_pipeline,
                  out_dir_fold = out_dir_fold)
  }

  lgr::lgr$info("===== 06 AGGREGATING RESULTS ACROSS FOLDS =====")
  aggregate_results(config = config, out_dir = out_dir)
  
  if (isTRUE(config$ensemble_models_params$use)) {
  lgr::lgr$info("===== 07 CREATING STACKED ENSEMBLE =====")

  ensemble_models(config = config,
                   data_raw = data,
                   out_dir = out_dir)
  }
  
  lgr::lgr$info("===== RUN COMPLETE =====")
  return(TRUE)
}

main_pipeline <- function(config, training_data, testing_data, out_dir_fold) {
  
  lgr::lgr$info(paste0("Starting main pipeline for fold ", gsub("fold_", "", basename(out_dir_fold))))
  
  rank_probes_dir <- file.path(out_dir_fold, "02_rank_probes")
  select_probes_dir <- file.path(out_dir_fold, "03_select_probes")
  train_models_dir <- file.path(out_dir_fold, "04_train_models")
  
  lgr::lgr$info("===== 02 RANKING PROBES =====")
  rank_probes(config = config, data = training_data, out_dir = out_dir_fold)

  lgr::lgr$info("===== 03 SELECTING PROBES =====")
  select_probes(config = config, data = training_data, rank_probes_dir = rank_probes_dir, out_dir = out_dir_fold)

  lgr::lgr$info("===== 04 TRAINING MODELS =====")
  train_models(config = config, data = training_data, filter_probes_dir = select_probes_dir, out_dir = out_dir_fold)

  lgr::lgr$info("===== 05 TESTING MODELS =====")
  test_models(config = config, data = testing_data, train_models_dir = train_models_dir, out_dir = out_dir_fold)
  
  lgr::lgr$info(paste0("Main pipeline complete for fold ", gsub("fold_", "", basename(out_dir_fold))))
  return(TRUE)
}