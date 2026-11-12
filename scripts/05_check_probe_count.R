# 05_check_probe_count.R
library(caret)
library(xgboost)

ens_path <- "./output/02_generate_mps/ensemble/ens_cortisol_auc.rds"
main_out_dir <- file.path("./output/05_check_probe_count")
ens <- readRDS(ens_path)

if (!dir.exists(main_out_dir)) dir.create(main_out_dir, recursive = TRUE)

sink(file.path(main_out_dir, "probe_counts.txt"))

model_features <- list()
all_coefnames <- list()

used_probes_df <- data.frame(Model = character(), Used_Probes = integer(), stringsAsFactors = FALSE)

for (model_name in names(ens$base_models)) {
  base_model <- ens$base_models[[model_name]]$model
  
  # Get base model probe importance
  importance_dt <- xgboost::xgb.importance(model = base_model$finalModel)
  used_features <- importance_dt$Feature
  
  model_features[[model_name]] <- used_features
  all_coefnames[[model_name]] <- base_model$coefnames
  
  # Get the number of probes available during model training
  total_passed_in <- length(base_model$coefnames)
  
  # Add the model name and number of used probes to df
  used_probes_df <- rbind(used_probes_df, data.frame(
    Model = model_name, 
    Used_Probes = length(used_features)
  ))
  
  # Print individual model statistics
  cat(sprintf("Model: %s\n", model_name))
  cat(sprintf("  Probes available during training: %d\n", total_passed_in))
  cat(sprintf("  Probes used by the model: %d\n", length(used_features)))
  cat(sprintf("  Probes not used by the model: %d\n\n", total_passed_in - length(used_features)))
}

# Get probe overlap
all_used_features <- unique(unlist(model_features))
total_training_probes <- length(unique(unlist(all_coefnames)))

cat(sprintf("Unique probes available to base models during training: %d\n", total_training_probes))
cat(sprintf("Unique probes after xgboost feature selection: %d\n\n", length(all_used_features)))

sink()

write.csv(used_probes_df, file.path(main_out_dir, "used_probe_counts.csv"), row.names = FALSE)