# 03_generate_predictions.R
library(caret)
library(lgr)
library(EpiDISH)

source("./scripts/generate_mps/01_process_data.R")

ens_path <- file.path("./output/02_generate_mps/ensemble/ens_cortisol_auc.rds")
mcig_cpgs <- read.csv("./data/mcigarette/Predictors_long.csv", header = TRUE)
GSE77445_formatted_beta_path <- "./output/01_prepare_data/GSE77445_formatted_beta.rds"
champagne_pre_formatted_beta_path <- "./output/01_prepare_data/champagne_pre_formatted_beta.rds"
champagne_post_formatted_beta_path <- "./output/01_prepare_data/champagne_post_formatted_beta.rds"

imputation_ref_path <- "./output/01_prepare_data/GSE119078_formatted_beta.rds"

main_out_dir <- file.path("./output/03_generate_predictions")

# Manually remove samples corresponding with participants that were identified as methylation outliers in either pre or post
champagne_pre_samples_to_remove <- c("205851990031_R07C01", "205837550052_R03C01")
champagne_post_samples_to_remove <- c("205837550052_R04C01", "205851990031_R08C01")

if (!dir.exists(main_out_dir)) dir.create(main_out_dir, recursive = TRUE)

# Generate mCigarette scores
generate_mcig_scores <- function(data, cpgs) {
  message(paste0("Generating mCigarette scores for ", deparse(substitute(data))))
  
  .na_to_mean <-function(methyl){
    methyl[is.na(methyl)] <- mean(methyl, na.rm = TRUE)
    return(methyl)
  }
  
  methyl <- data$methyl
  pheno <- data$pheno
  sites <- c(cpgs$CpG)
  required <- c(cpgs[cpgs$Predictor == "mCigarette", "CpG"])
  
  coef <- methyl[intersect(rownames(methyl), sites), ]
  cols <- colnames(coef)
  coef <- as.data.frame(t(apply(coef, 1, scale)))
  colnames(coef) <- cols
  coef <- as.data.frame(t(apply(coef, 1, function(x) .na_to_mean(x))))
  
  tmp <- coef[intersect(row.names(coef), cpgs[cpgs$Predictor == "mCigarette", "CpG"]), ]
  tmp_coef <- cpgs[cpgs$Predictor == "mCigarette", ]
  tmp_coef <- tmp_coef[match(row.names(tmp), tmp_coef$CpG), ]
  
  out <- data.frame()
  out[colnames(coef), "mCigarette"] <- colSums(tmp_coef$Beta * tmp)
  
  out$ID <- row.names(out)
  pheno$ID <- row.names(pheno)
  pheno <- merge(pheno, out, by = "ID", all.y = TRUE, sort = FALSE)
  rownames(pheno) <- pheno$ID
  
  # Re-align methyl columns to match pheno rows
  common_ids <- intersect(rownames(pheno), colnames(data$methyl))
  pheno <- pheno[common_ids, , drop = FALSE]
  methyl <- data$methyl[, common_ids, drop = FALSE]
  
  list(pheno = pheno, methyl = methyl)
}
                                
# Get cell proportion estimates
data(centEpiFibIC.m)
data(centBloodSub.m)
data(cent12CT450k.m)

apply_hepidish <- function(data_obj) {
  message("Estimating cell proportions (hierarchical)...")
  prop <- EpiDISH::hepidish(beta.m = data_obj$methyl, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = "RPC")
  prop <- as.data.frame(prop)
  
  data_obj$pheno$Epi <- prop$Epi
  data_obj$pheno$Neutro <- prop$Neutro
  
  prop_leuko <- subset(prop, select = -c(Epi, Fib))
  data_obj$pheno$Leuko_sum <- rowSums(prop_leuko)
  
  prop_lympho <- subset(prop, select = c(B, NK, CD4T, CD8T))
  data_obj$pheno$lympho_sum <- rowSums(prop_lympho)
  
  cell_subset <- subset(prop, select = c(Neutro, Eosino))
  data_obj$pheno$Gran <- rowSums(cell_subset)

  data_obj$pheno$nlr <- data_obj$pheno$Neutro / (data_obj$pheno$lympho_sum + 1e-9)
    
  # Add specific immune cells for later covariate adjustment
  data_obj$pheno$B <- prop$B
  data_obj$pheno$NK <- prop$NK
  data_obj$pheno$CD4T <- prop$CD4T
  data_obj$pheno$CD8T <- prop$CD8T
  data_obj$pheno$Mono <- prop$Mono
  
  return(data_obj)
}

apply_epidish_blood <- function(data_obj) {
  message("Estimating cell proportions (blood)...")
  prop <- EpiDISH::epidish(beta.m = data_obj$methyl, ref.m = cent12CT450k.m, method = "RPC")
  prop <- as.data.frame(prop$estF)
  cell_subset <- subset(prop, select = c(Bmem, Bnv)); data_obj$pheno$B <- rowSums(cell_subset)
  cell_subset <- subset(prop, select = c(CD4Tmem, CD4Tnv, Treg)); data_obj$pheno$CD4T <- rowSums(cell_subset)
  cell_subset <- subset(prop, select = c(CD8Tmem, CD8Tnv)); data_obj$pheno$CD8T <- rowSums(cell_subset)
  cell_subset <- subset(prop, select = c(Neu, Baso, Eos)); data_obj$pheno$Gran <- rowSums(cell_subset)
  data_obj$pheno$NK <- prop$NK
  data_obj$pheno$Mono <- prop$Mono
  data_obj$pheno$Neutro <- prop$Neu
  prop_lympho <- subset(data_obj$pheno, select = c(B, NK, CD4T, CD8T))
  data_obj$pheno$lympho_sum <- rowSums(prop_lympho)
  data_obj$pheno$nlr <- data_obj$pheno$Neutro / (data_obj$pheno$lympho_sum + 1e-9)
  return(data_obj)
}

# Get predictions from previously trained model
generate_predictions <- function(new_data, model_path, out_dir, target_pheno = NULL, samples_to_remove = NULL, seed = 123, imputation_ref_path = NULL) {

  log_filename <- paste0("log.txt")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  log_filepath <- file.path(out_dir, log_filename)
  lgr::lgr$add_appender(lgr::AppenderFile$new(log_filepath), name = "file_log")
  on.exit({ lgr::lgr$remove_appender("file_log") }, add = TRUE)
  
  data_set_name <- deparse(substitute(new_data))
  lgr::lgr$info("Starting prediction generation...")
  
  if (!file.exists(model_path)) stop("Model not found at: ", model_path)
  model_object <- readRDS(model_path)
  model_name <- tools::file_path_sans_ext(basename(model_path))
  config <- model_object$metadata$config
  is_ensemble <- "meta_model" %in% names(model_object)
  
  has_phenotype <- !is.null(target_pheno)
  phenotype_name <- NULL
  if (has_phenotype) {
    config$process_data_params$pheno$pheno_name <- target_pheno
    phenotype_name <- config$process_data_params$pheno$pheno_name
  } else {
    if (!is.null(config$process_data_params$pheno)) config$process_data_params$pheno$process <- FALSE
  }
  
  required_probes <- if(is_ensemble) unique(unlist(lapply(model_object$base_models, function(m) m$model$coefnames))) else model_object$model$coefnames
  
  present_probes <- rownames(new_data$methyl)
  missing_probes <- setdiff(required_probes, present_probes)
  
  if (length(missing_probes) > 0) {
    lgr::lgr$info(sprintf("Missing %d required probes. Imputing using reference...", length(missing_probes)))
    
    if (is.null(imputation_ref_path) || !file.exists(imputation_ref_path)) {
      stop("Error: Imputation reference could not be found.")
    }
    
    ref_data <- readRDS(imputation_ref_path)
    if (!"methyl" %in% names(ref_data)) stop("Error: Imputation reference does not contain methyl.")
    
    # Identify which missing probes are in reference compared to those that are completely missing
    available_in_ref <- intersect(missing_probes, rownames(ref_data$methyl))
    still_missing <- setdiff(missing_probes, available_in_ref)
    
    imputed_csv_filename <- paste0("imputed_probes.csv")
    imputed_csv_path <- file.path(out_dir, imputed_csv_filename)
    
    imputed_summary <- data.frame(
      Probe_ID = c(available_in_ref, still_missing),
      Imputation_Method = c(rep("Reference_Median", length(available_in_ref)),
                            rep("Global_Fallback", length(still_missing)))
    )
    
    write.csv(imputed_summary, imputed_csv_path, row.names = FALSE)
    lgr::lgr$info(paste("Saved list of imputed probes to:", imputed_csv_path))
    
    n_samples <- ncol(new_data$methyl)
    imputed_chunks <- list()
    
    if (length(available_in_ref) > 0) {
      lgr::lgr$info(sprintf("Imputing %d probes using probe-wise reference medians...", length(available_in_ref)))
      ref_subset <- ref_data$methyl[available_in_ref, , drop = FALSE]
      
      ref_values <- apply(ref_subset, 1, median, na.rm = TRUE)
      
      chunk_matrix <- matrix(rep(ref_values, times = n_samples), 
                             nrow = length(ref_values), 
                             ncol = n_samples, 
                             byrow = FALSE)
      rownames(chunk_matrix) <- names(ref_values)
      imputed_chunks[[1]] <- chunk_matrix
    }
    
    if (length(still_missing) > 0) {
      lgr::lgr$info(sprintf("%d probes are missing from both new data and reference data.", length(still_missing)))
      lgr::lgr$info("Calculating global median fallback from all available required probes in reference...")
      
      global_pool_probes <- intersect(required_probes, rownames(ref_data$methyl))
      
      if (length(global_pool_probes) == 0) {
        stop("Error: Imputation reference contains no required probes.")
      }
      
      ref_global_subset <- ref_data$methyl[global_pool_probes, , drop = FALSE]
      global_fallback_value <- median(as.numeric(ref_global_subset), na.rm = TRUE)
      
      lgr::lgr$info(sprintf("Imputing %d probes that are completely missing from new data and the reference using global median of available and required probes in the reference: %f", length(still_missing), global_fallback_value))
      
      chunk_matrix_fallback <- matrix(global_fallback_value, 
                                      nrow = length(still_missing), 
                                      ncol = n_samples)
      rownames(chunk_matrix_fallback) <- still_missing
      imputed_chunks[[2]] <- chunk_matrix_fallback
    }
    
    # Combine and append
    if (length(imputed_chunks) > 0) {
      full_imputed_matrix <- do.call(rbind, imputed_chunks)
      colnames(full_imputed_matrix) <- colnames(new_data$methyl)
      
      new_data$methyl <- rbind(new_data$methyl, full_imputed_matrix)
      lgr::lgr$info("Imputation complete.")
    }
    
  } else {
    lgr::lgr$info("All required probes present. No imputation needed.")
  }
  
  results_dir <- file.path(out_dir, data_set_name)
  if (!dir.exists(results_dir)) dir.create(results_dir)
  
  process_data_result <- process_data(config = config, data = new_data, out_dir = results_dir, process_pheno = has_phenotype, required_probes = required_probes, samples_to_remove = samples_to_remove)
  processed_data <- readRDS(file.path(process_data_result$dir, "processed_data.rds"))
  
  present_probes_proc <- intersect(required_probes, rownames(processed_data$methyl))
  if(length(present_probes_proc) != length(required_probes)) {
    stop("Error: Discrepancy in probe counts after processing, despite imputation.")
  }
  
  final_predictions <- NULL
  base_predictions_df <- NULL
  
  if (is_ensemble) {
    base_models <- model_object$base_models
    
    base_predictions_df <- as.data.frame(lapply(seq_along(base_models), function(i) {
      bm <- base_models[[i]]$model
      bm_name <- names(base_models)[i]
      if(is.null(bm_name)) bm_name <- paste0("Model_", i)
      
      bm_probes <- bm$coefnames
      
      pred_df <- as.data.frame(t(processed_data$methyl[bm_probes, , drop = FALSE]))
      predict(bm, newdata = pred_df)
    }))
    
    if (!is.null(names(base_models))) colnames(base_predictions_df) <- names(base_models)
    final_predictions <- predict(model_object$meta_model, newdata = base_predictions_df)    
  } else {
    stop("Error: Model should be ensemble.")
  }
  
  final_samples <- colnames(processed_data$methyl)
  names(final_predictions) <- final_samples
  
  output_df <- if(has_phenotype) {
    df <- new_data$pheno[final_samples, , drop = FALSE]
    df[[paste0(phenotype_name, "_proc")]] <- processed_data$pheno[final_samples, phenotype_name]
    df[[paste0(phenotype_name, "_pred")]] <- final_predictions
    df
  } else {
    df <- if(!is.null(new_data$pheno) && nrow(new_data$pheno) > 0) new_data$pheno[final_samples, , drop = FALSE] else data.frame(row.names = final_samples)
    df$pred <- final_predictions
    df
  }
  
  if (is_ensemble && !is.null(base_predictions_df)) {
    rownames(base_predictions_df) <- final_samples
    output_df <- cbind(output_df, base_predictions_df)
  }
  
  processed_data$pheno <- output_df
  rds_output_filename <- file.path(results_dir, paste0(model_name, "_on_", data_set_name, "_proc.rds"))
  saveRDS(processed_data, rds_output_filename)
  
  file.remove(file.path(process_data_result$dir, "processed_data.rds"))
  return(processed_data)
}
                                                          
message("Getting predictions...")
champagne_pre_formatted_beta <- readRDS(champagne_pre_formatted_beta_path)
champagne_post_formatted_beta <- readRDS(champagne_post_formatted_beta_path)
GSE77445_formatted_beta <- readRDS(GSE77445_formatted_beta_path)
                                                          
champagne_pre_pred <- generate_mcig_scores(data = champagne_pre_formatted_beta, cpgs = mcig_cpgs)
champagne_post_pred <- generate_mcig_scores(data = champagne_post_formatted_beta, cpgs = mcig_cpgs)
GSE77445_pred <- generate_mcig_scores(data = GSE77445_formatted_beta, cpgs = mcig_cpgs)
                                                          
champagne_pre_pred <- apply_hepidish(champagne_pre_pred)
champagne_post_pred <- apply_hepidish(champagne_post_pred)
GSE77445_pred <- apply_epidish_blood(GSE77445_pred)
  
champagne_pre_pred <- generate_predictions(new_data = champagne_pre_pred, model_path = ens_path, target_pheno = "AUCi_scc", out_dir = main_out_dir, samples_to_remove = champagne_pre_samples_to_remove, imputation_ref_path = imputation_ref_path)
champagne_post_pred <- generate_predictions(new_data = champagne_post_pred, model_path = ens_path, target_pheno = "AUCi_scc", out_dir = main_out_dir, samples_to_remove = champagne_post_samples_to_remove, imputation_ref_path = imputation_ref_path)
GSE77445_pred <- generate_predictions(new_data = GSE77445_pred, model_path = ens_path, target_pheno = "cortisol_auc", out_dir = main_out_dir, samples_to_remove = NULL, imputation_ref_path = imputation_ref_path)
                                                          
write.csv(champagne_pre_pred$pheno, file.path(main_out_dir, "champagne_pre_pred.csv"), row.names = FALSE)
write.csv(champagne_post_pred$pheno, file.path(main_out_dir, "champagne_post_pred.csv"), row.names = FALSE)
write.csv(GSE77445_pred$pheno, file.path(main_out_dir, "GSE77445_pred.csv"), row.names = FALSE)