# 07_generate_simple_models.R
library(lgr)
library(caret)
library(glmnet)
library(dplyr)
library(tidyr)
library(boot)
library(purrr)

source("./scripts/generate_mps/utils.R")
source("./scripts/generate_mps/01_process_data.R")

champagne_pre_data_path <- file.path("./output/01_prepare_data/champagne_pre_formatted_beta.rds")
champagne_post_data_path <- file.path("./output/01_prepare_data/champagne_post_formatted_beta.rds")
champagne_pre_pred_path <- file.path("./output/06_get_adjusted_predictions/champagne_pre_pred_adj.csv")
champagne_post_pred_path <- file.path("./output/06_get_adjusted_predictions/champagne_post_pred_adj.csv")
config <- readRDS("./data/config.rds")
main_out_dir <- file.path("./output/07_generate_simple_models")

if (!dir.exists(main_out_dir)) dir.create(main_out_dir, recursive = TRUE)

log_filepath <- file.path(main_out_dir, paste0("log.txt"))
lgr::lgr$add_appender(lgr::AppenderFile$new(log_filepath), name = "file_log")

set.seed(config$main_pipeline_params$seed)

lgr::lgr$info("Processing training data...")
gse77445_raw <- readRDS("./output/01_prepare_data/GSE77445_formatted_beta.rds")

train_proc_result <- process_data(
  config = config, 
  data = gse77445_raw, 
  out_dir = main_out_dir, 
  process_methyl = TRUE, 
  process_pheno = TRUE, 
  sub_dir_name = "train_processed"
)

train_proc_data <- readRDS(file.path(train_proc_result$dir, "processed_data.rds"))
train_methyl <- t(train_proc_data$methyl)
train_pheno <- train_proc_data$pheno[[config$process_data_params$pheno$pheno_name]]

lgr::lgr$info("Selecting top 850 probes by absolute correlation...")
cors <- cor(train_methyl, train_pheno, use = "pairwise.complete.obs")
abs_cors <- abs(cors[, 1])
top_850_probes <- names(sort(abs_cors, decreasing = TRUE)[1:850])
train_methyl_850 <- train_methyl[, top_850_probes, drop = FALSE]

lgr::lgr$info("Filtering collinear probes...")
cor_mat_850 <- cor(train_methyl_850, use = "pairwise.complete.obs")
cor_mat_850[is.na(cor_mat_850)] <- 0 
high_cor_idx <- caret::findCorrelation(cor_mat_850, cutoff = 0.9, exact = FALSE)

if (length(high_cor_idx) > 0) {
  train_methyl_filtered <- train_methyl_850[, -high_cor_idx, drop = FALSE]
} else {
  train_methyl_filtered <- train_methyl_850
}

final_train_probes <- colnames(train_methyl_filtered)
lgr::lgr$info(sprintf("Retained %d probes after correlation filtering.", length(final_train_probes)))

lgr::lgr$info("Training lasso and elastic net models...")

train_ctrl <- caret::trainControl(
  method = "repeatedcv", 
  number = 5,
  repeats = 10,
  summaryFunction = r_summary,
  allowParallel = FALSE
)

lasso_grid <- expand.grid(alpha = 1, lambda = seq(0.0001, 1, length.out = 100))
lasso_model <- caret::train(
  x = train_methyl_filtered, 
  y = train_pheno,
  method = "glmnet",
  metric = "Rsquared",
  tuneGrid = lasso_grid,
  trControl = train_ctrl
)

enet_grid <- expand.grid(alpha = seq(0.1, 0.9, by = 0.1), lambda = seq(0.0001, 1, length.out = 100))
enet_model <- caret::train(
  x = train_methyl_filtered, 
  y = train_pheno,
  method = "glmnet",
  metric = "Rsquared",
  tuneGrid = enet_grid,
  trControl = train_ctrl
)

saveRDS(lasso_model, file.path(main_out_dir, "simple_lasso_model.rds"))
saveRDS(enet_model, file.path(main_out_dir, "simple_enet_model.rds"))

required_probes <- final_train_probes

lgr::lgr$info("Processing test cohorts and generating predictions...")

imputation_ref <- readRDS("./output/01_prepare_data/GSE119078_formatted_beta.rds")

# Simple prediction function for the lasso and elastic net models
generate_simple_predictions <- function(raw_data_path, dataset_name, required_probes) {
  test_data <- readRDS(raw_data_path)
  
  missing_probes <- setdiff(required_probes, rownames(test_data$methyl))
  if (length(missing_probes) > 0) {
    lgr::lgr$info(sprintf("Imputing %d missing probes for %s...", length(missing_probes), dataset_name))
    
    available_in_ref <- intersect(missing_probes, rownames(imputation_ref$methyl))
    still_missing <- setdiff(missing_probes, available_in_ref)
    n_samples <- ncol(test_data$methyl)
    imputed_chunks <- list()
    
    if (length(available_in_ref) > 0) {
      ref_subset <- imputation_ref$methyl[available_in_ref, , drop = FALSE]
      ref_values <- apply(ref_subset, 1, median, na.rm = TRUE)
      chunk_matrix <- matrix(rep(ref_values, times = n_samples), nrow = length(ref_values), ncol = n_samples, byrow = FALSE)
      rownames(chunk_matrix) <- names(ref_values)
      imputed_chunks[[1]] <- chunk_matrix
    }
    
    if (length(still_missing) > 0) {
      global_pool_probes <- intersect(required_probes, rownames(imputation_ref$methyl))
      ref_global_subset <- imputation_ref$methyl[global_pool_probes, , drop = FALSE]
      global_fallback_value <- median(as.numeric(ref_global_subset), na.rm = TRUE)
      chunk_matrix_fallback <- matrix(global_fallback_value, nrow = length(still_missing), ncol = n_samples)
      rownames(chunk_matrix_fallback) <- still_missing
      imputed_chunks[[2]] <- chunk_matrix_fallback
    }
    
    if (length(imputed_chunks) > 0) {
      full_imputed_matrix <- do.call(rbind, imputed_chunks)
      colnames(full_imputed_matrix) <- colnames(test_data$methyl)
      test_data$methyl <- rbind(test_data$methyl, full_imputed_matrix)
    }
  }
  
  # Process data using config
  proc_result <- process_data(
    config = config, 
    data = test_data, 
    out_dir = main_out_dir, 
    process_methyl = TRUE, 
    process_pheno = FALSE, 
    required_probes = required_probes, 
    sub_dir_name = paste0(dataset_name, "_processed")
  )
  proc_data <- readRDS(file.path(proc_result$dir, "processed_data.rds"))  
  pred_df <- as.data.frame(t(proc_data$methyl[required_probes, , drop = FALSE]))
  
  list(
    lasso = predict(lasso_model, newdata = pred_df),
    enet = predict(enet_model, newdata = pred_df)
  )
}

# Generate predictions
preds_pre <- generate_simple_predictions(champagne_pre_data_path, "champagne_pre", required_probes)
preds_post <- generate_simple_predictions(champagne_post_data_path, "champagne_post", required_probes)

# Get previously-generated predictions
champagne_pre_pred <- read.csv(champagne_pre_pred_path)
champagne_post_pred <- read.csv(champagne_post_pred_path)

# Add lasso and elastic net predictions
champagne_pre_pred$lasso_pred <- preds_pre$lasso[champagne_pre_pred$ID]
champagne_pre_pred$enet_pred <- preds_pre$enet[champagne_pre_pred$ID]

champagne_post_pred$lasso_pred <- preds_post$lasso[champagne_post_pred$ID]
champagne_post_pred$enet_pred <- preds_post$enet[champagne_post_pred$ID]

# Save the updated data
write.csv(champagne_pre_pred, file.path(main_out_dir, "champagne_pre_pred_glmnet.csv"), row.names = FALSE)
write.csv(champagne_post_pred, file.path(main_out_dir, "champagne_post_pred_glmnet.csv"), row.names = FALSE)