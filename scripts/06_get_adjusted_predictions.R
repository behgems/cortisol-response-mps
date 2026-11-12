# 06_get_adjusted_predictions
library(dplyr)
library(stats)

pre_proc_path <- "./output/03_generate_predictions/champagne_pre_pred/ens_cortisol_auc_on_champagne_pre_pred_proc.rds"
post_proc_path <- "./output/03_generate_predictions/champagne_post_pred/ens_cortisol_auc_on_champagne_post_pred_proc.rds"
imputed_csv_path <- "./output/03_generate_predictions/imputed_probes.csv"
aggregated_shap_path <- "./output/04_analyse_shap_importance/final_shap_results_aggregated.rds"
main_out_dir <- file.path("./output/06_get_adjusted_predictions")

dir.create(main_out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading data...")
champagne_pre_proc <- readRDS(pre_proc_path)
champagne_post_proc <- readRDS(post_proc_path)
imputed_probes <- read.csv(imputed_csv_path)$Probe_ID

# Get the SHAP results
final_shap_results <- readRDS(aggregated_shap_path)
shap_pre <- final_shap_results[["Champagne_Pre"]]
shap_post <- final_shap_results[["Champagne_Post"]]

# Get predictions adjusted for SHAP values of imputed probes
adjust_predictions <- function(dataset_name, shap_matrix, proc_obj) {
  pheno <- proc_obj$pheno
  shap_matrix <- shap_matrix[rownames(pheno), ] # Ensure matching order
  
  # Identify imputed probes
  valid_imputed <- intersect(imputed_probes, colnames(shap_matrix))
  
  # For each participant, sum up the SHAP values of the imputed probes
  imputed_shap_sum <- rowSums(shap_matrix[, valid_imputed, drop = FALSE])
  
  # Ensure ID column is present
  if (!"ID" %in% colnames(pheno)) {
    pheno$ID <- rownames(pheno)
  }
  
  # Add to pheno
  pheno$observed_proc <- pheno$AUCi_scc_proc
  pheno$imputed_shift  <- imputed_shap_sum
  pheno$AUCi_scc_pred_adj <- pheno$AUCi_scc_pred - imputed_shap_sum
  
  return(invisible(pheno))
}

message("Adjusting predictions...")
adj_pre <- adjust_predictions("champagne_pre", shap_pre, champagne_pre_proc)
adj_post <- adjust_predictions("champagne_post", shap_post, champagne_post_proc)

write.csv(adj_pre, file.path(main_out_dir, "champagne_pre_pred_adj.csv"), row.names = FALSE)
write.csv(adj_post, file.path(main_out_dir, "champagne_post_pred_adj.csv"), row.names = FALSE)