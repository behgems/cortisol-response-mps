# 08_run_primary_analyses.R
library(ggplot2)
library(tidyr)
library(dplyr)
library(fastDummies)
library(Hmisc)
library(lm.beta)
library(effsize)
library(viridis)
library(effectsize)
library(car)
library(boot)
library(purrr)

champagne_pre_pred <- read.csv("./output/07_generate_simple_models/champagne_pre_pred_glmnet.csv")
champagne_post_pred <- read.csv("./output/07_generate_simple_models/champagne_post_pred_glmnet.csv")
GSE77445_pred <- read.csv("./output/03_generate_predictions/GSE77445_pred.csv")
used_probes_df <- read.csv("./output/05_check_probe_count/used_probe_counts.csv")
main_out_dir <- file.path("./output/08_run_primary_analyses")

if (!dir.exists(main_out_dir)) dir.create(main_out_dir, recursive = TRUE)

# Set up factors
# Merge Other into Hispanic because there are too few participants in the former (n = 3)
champagne_pre_pred$Ethnicity_Race[champagne_pre_pred$Ethnicity_Race == "Other"] <- "Hispanic"
champagne_post_pred$Ethnicity_Race[champagne_post_pred$Ethnicity_Race == "Other"] <- "Hispanic"

champagne_pre_pred$Ethnicity_Race <- as.factor(champagne_pre_pred$Ethnicity_Race)
champagne_pre_pred$Ethnicity_Race <- relevel(champagne_pre_pred$Ethnicity_Race, ref = "Asian")
champagne_pre_pred <- dummy_cols(champagne_pre_pred, select_columns = "Ethnicity_Race", remove_first_dummy = FALSE)
champagne_pre_pred$Sex <- as.factor(champagne_pre_pred$Sex)
champagne_pre_pred$Sex <- relevel(champagne_pre_pred$Sex, ref = "female")
champagne_pre_pred <- dummy_cols(champagne_pre_pred, select_columns = "Sex", remove_first_dummy = TRUE)

champagne_post_pred$Ethnicity_Race <- as.factor(champagne_post_pred$Ethnicity_Race)
champagne_post_pred$Ethnicity_Race <- relevel(champagne_post_pred$Ethnicity_Race, ref = "Asian")
champagne_post_pred <- dummy_cols(champagne_post_pred, select_columns = "Ethnicity_Race", remove_first_dummy = FALSE)
champagne_post_pred$Sex <- as.factor(champagne_post_pred$Sex)
champagne_post_pred$Sex <- relevel(champagne_post_pred$Sex, ref = "female")
champagne_post_pred <- dummy_cols(champagne_post_pred, select_columns = "Sex", remove_first_dummy = TRUE)

GSE77445_pred$gender <- as.factor(GSE77445_pred$gender)
GSE77445_pred$gender <- relevel(GSE77445_pred$gender, ref = "Female")
GSE77445_pred <- dummy_cols(GSE77445_pred, select_columns = "gender", remove_first_dummy = TRUE)

# Assign sample IDs
rownames(champagne_pre_pred) <- champagne_pre_pred$ID
rownames(champagne_post_pred) <- champagne_post_pred$ID
rownames(GSE77445_pred) <- GSE77445_pred$ID

# Get correlation matrices
champagne_cols_to_select <- c("Age", "Epi", "B", "NK", "CD4T", "Mono", "Gran", "lympho_sum", "Neutro", "Leuko_sum", "nlr", "AUCi_scc_proc", "Sex_male", "mCigarette", "AUCi_scc_pred", "lasso_pred", "enet_pred", "Ethnicity_Race_Asian", "Ethnicity_Race_Hispanic", "s1")

GSE77445_cols_to_select <- c("cortisol_auc_proc", "age", "gender_Male", "mCigarette", "B", "NK", "CD4T", "Mono", "Gran", "lympho_sum", "Neutro", "nlr", "cortisol_auc_pred")

create_final_corr_matrix <- function(r_matrix, p_matrix) {
  format_p_vec <- function(p_vals) {
    formatted <- ifelse(p_vals < 0.001, 
                        "<0.001", 
                        round(p_vals, 3))
    return(formatted)
  }
  formatted_matrix <- matrix(
    data = NA_character_, 
    nrow = nrow(r_matrix), 
    ncol = ncol(r_matrix),
    dimnames = dimnames(r_matrix)
  )
  for (i in seq_len(ncol(formatted_matrix))) {
    r_vals <- r_matrix[, i]
    p_vals <- p_matrix[, i]
    formatted_matrix[, i] <- paste0(round(r_vals, 3), " (", format_p_vec(p_vals), ")")
  }
  final_df <- as.data.frame(formatted_matrix, check.names = FALSE)
  return(final_df)
}

subset_df <- champagne_pre_pred[, champagne_cols_to_select]
corr_matrix <- rcorr(as.matrix(subset_df), type = "pearson")
final_corr_matrix <- create_final_corr_matrix(r_matrix = corr_matrix$r, p_matrix = corr_matrix$P)
write.csv(final_corr_matrix, file.path(main_out_dir, "champagne_pre_final_corr.csv"), row.names = TRUE)

subset_df <- champagne_post_pred[, champagne_cols_to_select]
corr_matrix <- rcorr(as.matrix(subset_df), type = "pearson")
final_corr_matrix <- create_final_corr_matrix(r_matrix = corr_matrix$r, p_matrix = corr_matrix$P)
write.csv(final_corr_matrix, file.path(main_out_dir, "champagne_post_final_corr.csv"), row.names = TRUE)

subset_df <- GSE77445_pred[, GSE77445_cols_to_select]
corr_matrix <- rcorr(as.matrix(subset_df), type = "pearson")
final_corr_matrix <- create_final_corr_matrix(r_matrix = corr_matrix$r, p_matrix = corr_matrix$P)
write.csv(final_corr_matrix, file.path(main_out_dir, "GSE77445_final_corr.csv"), row.names = TRUE)

sink(file.path(main_out_dir, "descriptives.txt"), append = FALSE)
champagne_cols_to_select_desc <- c("Age", "Epi", "B", "NK", "CD4T", "CD8T", "Mono", "Gran", "Leuko_sum", "AUCi_scc_proc", "Sex", "mCigarette", "AUCi_scc_pred", "lasso_pred", "enet_pred", "Ethnicity_Race", "AUCi_scc", "s1")

GSE77445_cols_to_select <- c("cortisol_auc_proc", "age", "gender", "mCigarette", "B", "NK", "CD4T", "CD8T", "Mono", "Gran", "cortisol_auc_pred", "cortisol_auc")

get_descriptives <- function(data, var_names) {
  for (var_name in var_names) {
    var <- data[[var_name]]            
    if (is.factor(var)) {      
      cat("\n========================================================================\n")
      cat(" ", var_name, "\n")
      cat("========================================================================\n")
      
      cat("\n----- Counts -----\n")
      tbl <- table(var, useNA = "always")
      print(tbl)
      
      cat("\n----- Percentages -----\n")
      prop_tbl_col <- round(prop.table(tbl) * 100, 2)
      print(prop_tbl_col)      
    } else { # Numeric
      cat("\n========================================================================\n")
      cat(" ", var_name, "\n")
      cat("========================================================================\n")
      
      var <- as.numeric(var)
      
      summary_overall <- summary(var)
      print(summary_overall)
      sd_val <- round(sd(var, na.rm = TRUE), 2)
      cat("SD =", sd_val, "\n")      
    }
  }
}

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
cat(" champagne_pre \n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

get_descriptives(data = champagne_pre_pred, var_names = champagne_cols_to_select_desc)
                
cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
cat(" champagne_post \n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

get_descriptives(data = champagne_post_pred, var_names = champagne_cols_to_select_desc)

cat("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")
cat(" GSE77445 \n")
cat("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

get_descriptives(data = GSE77445_pred, var_names = GSE77445_cols_to_select)

sink()

# Check change from pre to post
sink(file.path(main_out_dir, "change_results.txt"), append = FALSE)
print("Change from pre to post stressor")
model <- t.test(champagne_pre_pred$AUCi_scc_pred, champagne_post_pred$AUCi_scc_pred, paired = TRUE)
print(model)
es <- cohen.d(champagne_post_pred$AUCi_scc_pred, champagne_pre_pred$AUCi_scc_pred, paired = TRUE)
print(es)
sink()

# Get base model columns
base_model_cols <- grep("^f_", colnames(champagne_pre_pred), value = TRUE)

if (length(base_model_cols) == 0) {
  stop("Could not find any base model columns starting with 'f_' in the data.")
}

# Define unadjusted models
model_pre_ens_unadjusted <- lm(AUCi_scc_proc ~ AUCi_scc_pred, data = champagne_pre_pred)
model_pre_ens_unadjusted_s1 <- lm(AUCi_scc_proc ~ AUCi_scc_pred + s1, data = champagne_pre_pred)
model_pre_ens_unadjusted_adj <- lm(AUCi_scc_proc ~ AUCi_scc_pred_adj, data = champagne_pre_pred)
model_pre_ens_unadjusted_raw <- lm(AUCi_scc ~ AUCi_scc_pred, data = champagne_pre_pred)

pre_base_models <- lapply(base_model_cols, function(col) {
  lm(as.formula(paste("AUCi_scc_proc ~", col)), data = champagne_pre_pred)
})
names(pre_base_models) <- paste0("model_pre_base", seq_along(base_model_cols), "_unadjusted")

model_post_ens_unadjusted <- lm(AUCi_scc_proc ~ AUCi_scc_pred, data = champagne_post_pred)
model_post_ens_unadjusted_s1 <- lm(AUCi_scc_proc ~ AUCi_scc_pred + s1, data = champagne_post_pred)
model_post_ens_unadjusted_adj <- lm(AUCi_scc_proc ~ AUCi_scc_pred_adj, data = champagne_post_pred)
model_post_ens_unadjusted_raw <- lm(AUCi_scc ~ AUCi_scc_pred, data = champagne_post_pred)

post_base_models <- lapply(base_model_cols, function(col) {
  lm(as.formula(paste("AUCi_scc_proc ~", col)), data = champagne_post_pred)
})
names(post_base_models) <- paste0("model_post_base", seq_along(base_model_cols), "_unadjusted")

# Define adjusted models
data <- champagne_pre_pred

# CD8T omitted because it was uniformly 0
model_pre_ens_fully_adjusted <- lm(
  data = data,
  AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran + AUCi_scc_pred
)

model_pre_ens_fully_adjusted_s1 <- lm(
  data = data,
  AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran + AUCi_scc_pred + s1
)

model_pre_ens_adjusted <- lm(
  data = data,
  AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + AUCi_scc_pred
)

model_pre_ens_adjusted_s1 <- lm(
  data = data,
  AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + AUCi_scc_pred + s1
)

model_pre_ens_adjusted_cell_only <- lm(
  data = data,
  AUCi_scc_proc ~ B + NK + CD4T + Mono + Gran + AUCi_scc_pred
)

model_pre_ens_adjusted_ethnicity_only <- lm(
  data = data,
  AUCi_scc_proc ~ Ethnicity_Race + AUCi_scc_pred
)

# Covariate-only
model_pre_covs_fully <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran)
model_pre_covs_demo <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette)

# Lasso and elastic net
model_pre_lasso_unadjusted <- lm(AUCi_scc_proc ~ lasso_pred, data = data)
model_pre_enet_unadjusted <- lm(AUCi_scc_proc ~ enet_pred, data = data)

# Fully adjusted
model_pre_lasso_fully_adjusted <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran + lasso_pred)
model_pre_enet_fully_adjusted <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran + enet_pred)

# Demographic and lifestyle only
model_pre_lasso_adjusted <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + lasso_pred)
model_pre_enet_adjusted <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + enet_pred)

# Simpler model plus ensemble
model_pre_ens_lasso_unadjusted <- lm(data = data, AUCi_scc_proc ~ AUCi_scc_pred + lasso_pred)
model_pre_ens_enet_unadjusted <- lm(data = data, AUCi_scc_proc ~ AUCi_scc_pred + enet_pred)

data <- champagne_post_pred

model_post_ens_fully_adjusted <- lm(
  data = data,
  AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran + AUCi_scc_pred
)

model_post_ens_fully_adjusted_s1 <- lm(
  data = data,
  AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran + AUCi_scc_pred + s1
)

model_post_ens_adjusted <- lm(
  data = data,
  AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + AUCi_scc_pred
)

model_post_ens_adjusted_s1 <- lm(
  data = data,
  AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + AUCi_scc_pred + s1
)

model_post_ens_adjusted_cell_only <- lm(
  data = data,
  AUCi_scc_proc ~ B + NK + CD4T + Mono + Gran + AUCi_scc_pred
)

model_post_ens_adjusted_ethnicity_only <- lm(
  data = data,
  AUCi_scc_proc ~ Ethnicity_Race + AUCi_scc_pred
)

model_post_covs_fully <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran)
model_post_covs_demo <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette)

model_post_lasso_unadjusted <- lm(AUCi_scc_proc ~ lasso_pred, data = data)
model_post_enet_unadjusted <- lm(AUCi_scc_proc ~ enet_pred, data = data)

model_post_lasso_fully_adjusted <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran + lasso_pred)
model_post_enet_fully_adjusted <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + B + NK + CD4T + Mono + Gran + enet_pred)

model_post_lasso_adjusted <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + lasso_pred)
model_post_enet_adjusted <- lm(data = data, AUCi_scc_proc ~ Sex + Age + Ethnicity_Race + mCigarette + enet_pred)

model_post_ens_lasso_unadjusted <- lm(data = data, AUCi_scc_proc ~ AUCi_scc_pred + lasso_pred)
model_post_ens_enet_unadjusted <- lm(data = data, AUCi_scc_proc ~ AUCi_scc_pred + enet_pred)

models <- c(
  list(
    model_pre_ens_unadjusted = model_pre_ens_unadjusted,
    model_post_ens_unadjusted = model_post_ens_unadjusted
  ),
  pre_base_models,
  post_base_models,
  list(
    model_pre_ens_adjusted = model_pre_ens_adjusted,
    model_post_ens_adjusted = model_post_ens_adjusted,
    model_pre_ens_fully_adjusted = model_pre_ens_fully_adjusted,
    model_post_ens_fully_adjusted = model_post_ens_fully_adjusted,
    model_pre_ens_adjusted_ethnicity_only = model_pre_ens_adjusted_ethnicity_only,
    model_post_ens_adjusted_ethnicity_only = model_post_ens_adjusted_ethnicity_only,
    model_pre_ens_adjusted_cell_only = model_pre_ens_adjusted_cell_only,
    model_post_ens_adjusted_cell_only = model_post_ens_adjusted_cell_only,
	model_pre_ens_unadjusted_raw = model_pre_ens_unadjusted_raw,
	model_post_ens_unadjusted_raw = model_post_ens_unadjusted_raw,
    model_pre_ens_unadjusted_s1 = model_pre_ens_unadjusted_s1,
	model_post_ens_unadjusted_s1 = model_post_ens_unadjusted_s1,
	model_pre_ens_adjusted_s1 = model_pre_ens_adjusted_s1,
    model_post_ens_adjusted_s1 = model_post_ens_adjusted_s1,
    model_pre_ens_fully_adjusted_s1 = model_pre_ens_fully_adjusted_s1,
    model_post_ens_fully_adjusted_s1 = model_post_ens_fully_adjusted_s1,
    model_pre_lasso_unadjusted = model_pre_lasso_unadjusted,
    model_post_lasso_unadjusted = model_post_lasso_unadjusted,
    model_pre_enet_unadjusted = model_pre_enet_unadjusted,
    model_post_enet_unadjusted = model_post_enet_unadjusted,
    model_pre_lasso_adjusted = model_pre_lasso_adjusted,
    model_post_lasso_adjusted = model_post_lasso_adjusted,
    model_pre_enet_adjusted = model_pre_enet_adjusted,
    model_post_enet_adjusted = model_post_enet_adjusted,
    model_pre_lasso_fully_adjusted = model_pre_lasso_fully_adjusted,
    model_post_lasso_fully_adjusted = model_post_lasso_fully_adjusted,
    model_pre_enet_fully_adjusted = model_pre_enet_fully_adjusted,
    model_post_enet_fully_adjusted = model_post_enet_fully_adjusted,
    model_pre_ens_lasso_unadjusted = model_pre_ens_lasso_unadjusted,
    model_post_ens_lasso_unadjusted = model_post_ens_lasso_unadjusted,
    model_pre_ens_enet_unadjusted = model_pre_ens_enet_unadjusted,
    model_post_ens_enet_unadjusted = model_post_ens_enet_unadjusted,
    model_pre_ens_unadjusted_adj = model_pre_ens_unadjusted_adj,
    model_post_ens_unadjusted_adj = model_post_ens_unadjusted_adj
  )
)

model_names <- names(models)

get_boot_stats_std <- function(model, name, R = 10000) {
  
  # Add relevant data source to model object
  if (grepl("_pre_", name)) {
    data_source <- champagne_pre_pred
  } else if (grepl("_post_", name)) {
    data_source <- champagne_post_pred
  }
  model$call$data <- quote(data_source)
  
  message(paste("Bootstrapping standardised and unstandardised estimates for:", name))
  
  boot_std_manual <- function(data, indices) {     
    d <- data[indices, ]
    
    fit <- try(lm(formula(model), data = d), silent = TRUE) 
    
    if (inherits(fit, "try-error")) {
      # Return NA vector of length equal to coefficients + unstandardised coefficients + R squared + Spearman's rho
      return(rep(NA, (length(coef(model)) * 2) + 2)) 
    }
    
    # Get unstandardised coefficients
    b <- coef(fit)
    
    if (length(b) != length(coef(model))) {
      b_aligned <- rep(NA, length(coef(model)))
      names(b_aligned) <- names(coef(model))
      common <- intersect(names(b), names(b_aligned))
      b_aligned[common] <- b[common]
      b <- b_aligned
    }
    
    b_unstd <- b
    names(b_unstd) <- paste0(names(b), "_Unstd")
    
    # Get standardised coefficients
    mf <- model.frame(formula(model), data = d, na.action = na.pass)
    mm <- model.matrix(formula(model), data = mf)

    # Get SD
    sd_x <- apply(mm, 2, sd, na.rm = TRUE)
    
    y_name <- all.vars(formula(model))[1]
    sd_y <- sd(d[[y_name]], na.rm = TRUE)
    b_std <- b * (sd_x / sd_y)
    r_sq <- summary(fit)$r.squared
    
    # Get Spearman's rho for observed values with predicted values
    observed_y <- d[[y_name]]
    pred_y <- predict(fit, newdata = d)
    spearman_rho <- cor(observed_y, pred_y, method = "spearman", use = "complete.obs")
    
    return(c(b_std, b_unstd, "R_squared" = r_sq, "Spearmans_rho" = spearman_rho))
  }
  
  tryCatch({
    boot_obj <- boot::boot(data = data_source, statistic = boot_std_manual, R = R)
    
    boot_se <- apply(boot_obj$t, 2, sd, na.rm = TRUE) 
    
    # Get bootstrapped CIs
    get_ci <- function(idx) {
      vec <- boot_obj$t[, idx]
      sd_val <- sd(vec, na.rm = TRUE)
      if (is.na(sd_val) || sd_val == 0) {
        return(c(NA, NA))
      }       
      boot.ci(boot_obj, index = idx, type = "bca")$bca[4:5]
    }
    
    n_params <- length(boot_obj$t0)
    ci_list <- lapply(1:n_params, get_ci)
    ci_mat <- do.call(rbind, ci_list)
    
    # Summary table
    df_res <- data.frame(
      Model_Name = name,
      Term = names(boot_obj$t0),
      Estimate = as.numeric(boot_obj$t0), 
      Boot_SE = as.numeric(boot_se),
      Boot_Lower_95 = ci_mat[, 1],
      Boot_Upper_95 = ci_mat[, 2]
    )
    
    # Get significance based on CI
    df_res$Significant_CI <- (df_res$Boot_Lower_95 > 0 & df_res$Boot_Upper_95 > 0) | 
                             (df_res$Boot_Lower_95 < 0 & df_res$Boot_Upper_95 < 0)
    
    df_res$Significant_CI[is.na(df_res$Significant_CI)] <- FALSE
    
    return(df_res)
    
  }, error = function(e) {
    message(paste("Error in", name, ":", e$message))
    return(NULL)
  })
}

set.seed(1995)
boot_results_std <- purrr::map2_dfr(models, names(models), get_boot_stats_std)
write.csv(boot_results_std, file.path(main_out_dir, "main_results.csv"), row.names = FALSE)

# Plotting
create_faceted_plot_from_stats <- function(df, stats_df, model_map, observed_col, filename, used_probes_map) {
    
  predicted_cols <- names(model_map)
  
  pretty_labels <- sapply(seq_along(predicted_cols), function(i) {
    col <- predicted_cols[i]
    if (col == "AUCi_scc_pred") return("Ensemble")
    if (col == "lasso_pred") return("Lasso")
    if (col == "enet_pred") return("Elastic Net")
    
    parts <- strsplit(col, "_m_")[[1]]
    if (length(parts) == 2) {
      model_part <- parts[2]
      ranking_part <- sub("^f_[0-9]+_r_", "", parts[1])
      
      ranking_part <- sub("_[0-9]+$", "", ranking_part)
      
      ranking_clean <- gsub("mutual_information", "MI", ranking_part)
      ranking_clean <- gsub("lasso_stability_selection", "LSS", ranking_clean)
      ranking_clean <- gsub("_", " | ", ranking_clean)
      
      used_count <- used_probes_map[[col]]
      return(paste0("Base Model ", i - 1, " | ", model_part, " | ", ranking_clean, " | ", used_count))
    }
    return(paste("Base Model", i - 1, "|", col))
  })
  
  # Pivot raw data to long format
  df_long <- df %>%
    pivot_longer(
      cols = all_of(predicted_cols),
      names_to = "raw_predictor",
      values_to = "predicted_value"
    ) %>%
    mutate(
      facet_label = factor(raw_predictor, 
                           levels = predicted_cols, 
                           labels = pretty_labels)
    )
  
  # Prepare annotation and regression line data
  plot_metadata_list <- lapply(predicted_cols, function(col_name) {
    
    target_model_name <- model_map[[col_name]]
    
    # Get standardised stats for the Label
    stat_row_std <- stats_df %>%
      filter(Model_Name == target_model_name, Term == col_name)
    
    # Get unstandardised stats for the regression line
    slope_row <- stats_df %>%
      filter(Model_Name == target_model_name, Term == paste0(col_name, "_Unstd"))
    
    intercept_row <- stats_df %>%
      filter(Model_Name == target_model_name, Term == "(Intercept)_Unstd")
    
    if(nrow(stat_row_std) == 0) return(NULL)
    
    # Create label text
    label_text <- sprintf(
      "italic(r) == %.2f ~ '[95%% CI = %.2f, %.2f]'",
      stat_row_std$Estimate,
      stat_row_std$Boot_Lower_95,
      stat_row_std$Boot_Upper_95
    )
    
    idx <- which(predicted_cols == col_name)
    current_pretty_label <- pretty_labels[idx]
    
    data.frame(
      raw_predictor = col_name,
      stat_label = label_text,
      facet_label = factor(current_pretty_label, levels = pretty_labels),
      # Add regression line params
      slope = slope_row$Estimate,
      intercept = intercept_row$Estimate
    )
  })
  
  metadata_df <- do.call(rbind, plot_metadata_list)
  
  line_color <- viridis::viridis(1, option = "B", begin = 0)
  dot_color <- viridis::viridis(1, option = "B", begin = 0.75)
  strip_bg_color <- viridis::viridis(1, option = "B", begin = 0.125)
  
  p_faceted <- ggplot(df_long, aes(x = predicted_value, y = .data[[observed_col]])) +
    geom_point(color = dot_color) +  
    
    # Get regression line with the bootstrapped unstandardised coefficients
    geom_abline(data = metadata_df, 
                aes(slope = slope, intercept = intercept), 
                linewidth = 1, color = line_color) +
    
    geom_label(
      data = metadata_df,
      aes(label = stat_label),
      parse = TRUE,
      x = -Inf, y = -Inf,
      hjust = -0.04,  
      vjust = -0.2,   
      size = 3.5,
      linewidth = 0,
      fill = "white",
      alpha = 0,      
      inherit.aes = FALSE
    ) +
    
    facet_wrap(~ facet_label, scales = "free") +
    labs(
      title = NULL,
      x = "Predicted Cortisol Response",
      y = "Observed Cortisol Response"
    ) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 12, margin = margin(b = 12)),
      axis.title.y = element_text(size = 12, margin = margin(r = 12)),
      axis.title.x = element_text(size = 12, margin = margin(t = 12)),
      strip.background = element_rect(fill = strip_bg_color, color = strip_bg_color),
      strip.text = element_text(size = 10, face = "bold", color = "white"),
      legend.position = "none"
    )
  
  ggsave(filename, plot = p_faceted, width = 8, height = 7)
}

champagne_pre_map <- list()
champagne_pre_map[["AUCi_scc_pred"]] <- "model_pre_ens_unadjusted"

for (i in seq_along(base_model_cols)) {
  champagne_pre_map[[base_model_cols[i]]] <- paste0("model_pre_base", i, "_unadjusted")
}

champagne_post_map <- list()
champagne_post_map[["AUCi_scc_pred"]] <- "model_post_ens_unadjusted"

for (i in seq_along(base_model_cols)) {
  champagne_post_map[[base_model_cols[i]]] <- paste0("model_post_base", i, "_unadjusted")
}

used_probes_map <- setNames(used_probes_df$Used_Probes, used_probes_df$Model)

create_faceted_plot_from_stats(
  df = champagne_pre_pred,
  stats_df = boot_results_std,
  model_map = champagne_pre_map,
  observed_col = "AUCi_scc_proc",
  used_probes_map = used_probes_map,
  filename = file.path(main_out_dir, "faceted_champagne_pre.pdf")
)

create_faceted_plot_from_stats(
  df = champagne_post_pred,
  stats_df = boot_results_std,
  model_map = champagne_post_map,
  observed_col = "AUCi_scc_proc",
  used_probes_map = used_probes_map,
  filename = file.path(main_out_dir, "faceted_champagne_post.pdf")
)

# Get Cook's D
cooks_pre_1 <- cooks.distance(model_pre_ens_unadjusted)
cooks_post_1 <- cooks.distance(model_post_ens_unadjusted)
cooks_pre_2 <- cooks.distance(model_pre_ens_adjusted)
cooks_post_2 <- cooks.distance(model_post_ens_adjusted)
cooks_pre_3 <- cooks.distance(model_pre_ens_fully_adjusted)
cooks_post_3 <- cooks.distance(model_post_ens_fully_adjusted)

cooks_df <- bind_rows(
  data.frame(
    ID = names(cooks_pre_1), 
    Model = "Pre-stress (Unadjusted)", 
    Cooks_Distance = as.numeric(cooks_pre_1)
  ),
  data.frame(
    ID = names(cooks_post_1), 
    Model = "Post-stress (Unadjusted)", 
    Cooks_Distance = as.numeric(cooks_post_1)
  ),
  data.frame(
    ID = names(cooks_pre_2), 
    Model = "Pre-stress (Adjusted)", 
    Cooks_Distance = as.numeric(cooks_pre_2)
  ),
  data.frame(
    ID = names(cooks_post_2), 
    Model = "Post-stress (Adjusted)", 
    Cooks_Distance = as.numeric(cooks_post_2)
  ),
    data.frame(
    ID = names(cooks_pre_3), 
    Model = "Pre-stress (Fully Adjusted)", 
    Cooks_Distance = as.numeric(cooks_pre_3)
  ),
  data.frame(
    ID = names(cooks_post_3), 
    Model = "Post-stress (Fully Adjusted)", 
    Cooks_Distance = as.numeric(cooks_post_3)
  )
)

write.csv(cooks_df, file.path(main_out_dir, "cooks_distances_ensemble.csv"), row.names = FALSE)

# Get R squared change (with paired bootstrapped CIs) to compare models
calc_delta_r2_boot <- function(base_mod, full_mod, comp_name, timepoint, R = 10000, conf = 0.95) {
  
  data_source <- if (timepoint == "pre") champagne_pre_pred else champagne_post_pred
  
  f_base <- formula(base_mod)
  f_full <- formula(full_mod)
  
  # Use the same complete cases for both models so R2 values are comparable
  vars <- unique(c(all.vars(f_base), all.vars(f_full)))
  d_cc <- data_source[complete.cases(data_source[, vars]), vars, drop = FALSE]
  
  if (nrow(d_cc) != nobs(base_mod) || nrow(d_cc) != nobs(full_mod)) {
    warning(paste0(comp_name, " (", timepoint, "): models were originally fitted on different ",
                   "numbers of cases; using ", nrow(d_cc), " complete cases for both."))
  }
  
  message(paste("Bootstrapping R2 change for:", comp_name, "|", timepoint))
  
  delta_stat <- function(data, indices) {
    d <- data[indices, , drop = FALSE]
    fit_base <- tryCatch(lm(f_base, data = d), error = function(e) NULL)
    fit_full <- tryCatch(lm(f_full, data = d), error = function(e) NULL)
    if (is.null(fit_base) || is.null(fit_full)) return(c(NA, NA, NA))
    r2_base <- summary(fit_base)$r.squared
    r2_full <- summary(fit_full)$r.squared
    c(R2_Base = r2_base, R2_Full = r2_full, Delta_R2 = r2_full - r2_base)
  }
  
  boot_obj <- boot::boot(data = d_cc, statistic = delta_stat, R = R)
  
  n_failed <- sum(!complete.cases(boot_obj$t))
  if (n_failed > 0) warning(paste(comp_name, timepoint, ":", n_failed, "failed bootstrap resamples"))
  
  get_ci <- function(idx, type) {
    tryCatch({
      ci <- boot::boot.ci(boot_obj, conf = conf, index = idx, type = type)
      slot_name <- if (type == "perc") "percent" else type
      ci[[slot_name]][4:5]
    }, error = function(e) {
      warning(paste("CI failed for", comp_name, timepoint, "index", idx, ":", e$message))
      c(NA, NA)
    })
  }
  
  ci_base_bca  <- get_ci(1, "bca")
  ci_full_bca  <- get_ci(2, "bca")
  ci_delta_bca <- get_ci(3, "bca")
  
  data.frame(
    Timepoint = timepoint,
    Comparison = comp_name,
    N = nrow(d_cc),
    R2_Base = round(boot_obj$t0[1], 4),
    R2_Base_Lower_95 = round(ci_base_bca[1], 4),
    R2_Base_Upper_95 = round(ci_base_bca[2], 4),
    R2_Full = round(boot_obj$t0[2], 4),
    R2_Full_Lower_95 = round(ci_full_bca[1], 4),
    R2_Full_Upper_95 = round(ci_full_bca[2], 4),
    Delta_R2 = round(boot_obj$t0[3], 4),
    Delta_R2_Boot_SE = round(sd(boot_obj$t[, 3], na.rm = TRUE), 4),
    Delta_R2_Lower_95_BCa = round(ci_delta_bca[1], 4),
    Delta_R2_Upper_95_BCa = round(ci_delta_bca[2], 4),
    row.names = NULL
  )
}

delta_r2_results <- bind_rows(
  calc_delta_r2_boot(model_pre_covs_fully, model_pre_ens_fully_adjusted, "covs (inc. cells) vs. covs (inc. cells) + ensemble", "pre"),
  calc_delta_r2_boot(model_pre_covs_fully, model_pre_lasso_fully_adjusted, "covs (inc. cells) vs. covs (inc. cells) + lasso", "pre"),
  calc_delta_r2_boot(model_pre_covs_fully, model_pre_enet_fully_adjusted, "covs (inc. cells) vs. covs (inc. cells) + elastic net", "pre"),
  
  calc_delta_r2_boot(model_pre_covs_demo, model_pre_ens_adjusted, "covs vs. covs + ensemble", "pre"),
  calc_delta_r2_boot(model_pre_covs_demo, model_pre_lasso_adjusted, "covs vs. covs + lasso", "pre"),
  calc_delta_r2_boot(model_pre_covs_demo, model_pre_enet_adjusted, "covs vs. covs + elastic net", "pre"),
  
  calc_delta_r2_boot(model_pre_lasso_unadjusted, model_pre_ens_lasso_unadjusted, "lasso vs. lasso + ensemble", "pre"),
  calc_delta_r2_boot(model_pre_enet_unadjusted, model_pre_ens_enet_unadjusted, "elastic net vs. elastic net + ensemble", "pre"),
  
  calc_delta_r2_boot(model_post_covs_fully, model_post_ens_fully_adjusted, "covs (inc. cells) vs. covs (inc. cells) + ensemble", "post"),
  calc_delta_r2_boot(model_post_covs_fully, model_post_lasso_fully_adjusted, "covs (inc. cells) vs. covs (inc. cells) + lasso", "post"),
  calc_delta_r2_boot(model_post_covs_fully, model_post_enet_fully_adjusted, "covs (inc. cells) vs. covs (inc. cells) + elastic net", "post"),
  
  calc_delta_r2_boot(model_post_covs_demo, model_post_ens_adjusted, "covs vs. covs + ensemble", "post"),
  calc_delta_r2_boot(model_post_covs_demo, model_post_lasso_adjusted, "covs vs. covs + lasso", "post"),
  calc_delta_r2_boot(model_post_covs_demo, model_post_enet_adjusted, "covs vs. covs + elastic net", "post"),
  
  calc_delta_r2_boot(model_post_lasso_unadjusted, model_post_ens_lasso_unadjusted, "lasso vs. lasso + ensemble", "post"),
  calc_delta_r2_boot(model_post_enet_unadjusted, model_post_ens_enet_unadjusted, "elastic net vs. elastic net + ensemble", "post"),
  
  calc_delta_r2_boot(model_pre_ens_unadjusted, model_pre_ens_unadjusted_adj, "ensemble predicted values adjusted for imputed probes vs. original predicted values", "pre"),
  calc_delta_r2_boot(model_post_ens_unadjusted, model_post_ens_unadjusted_adj, "ensemble predicted values adjusted for imputed probes vs. original predicted values", "post")
)

write.csv(delta_r2_results, file.path(main_out_dir, "r_squared_change_results.csv"), row.names = FALSE)

# Compare GSE77445 with Champagne
# Champagne column to GSE77445 column
cohort_var_map <- list(
  Age           = "age",
  Sex           = "gender",
  mCigarette    = "mCigarette",
  B             = "B",
  NK            = "NK",
  CD4T          = "CD4T",
  CD8T          = "CD8T",
  Mono          = "Mono",
  Gran          = "Gran",
  AUCi_scc_proc = "cortisol_auc_proc",
  AUCi_scc_pred = "cortisol_auc_pred",
  AUCi_scc      = "cortisol_auc"
)

recode_sex_numeric <- function(x, label) {
  x_clean <- tolower(trimws(as.character(x)))
  bad_levels <- setdiff(unique(x_clean), c("female", "male"))
  if (length(bad_levels) > 0) {
    stop(paste0("Unexpected sex levels in ", label, ": ",
                paste(bad_levels, collapse = ", "),
                ". Expected only 'female' and 'male'."), call. = FALSE)
  }
  as.numeric(x_clean == "female")
}

get_cohort_values <- function(data, col, label) {
  if (!col %in% colnames(data)) {
    stop(paste0("Column '", col, "' not found in ", label, "."), call. = FALSE)
  }
  values <- data[[col]]
  values <- values[!is.na(values)]
  if (col %in% c("Sex", "gender")) {
    values <- recode_sex_numeric(values, paste0(label, "$", col))
  } else {
    values <- as.numeric(values)
  }
  if (length(values) < 3) {
    stop(paste0("Fewer than 3 non-missing values for '", col, "' in ", label, "."),
         call. = FALSE)
  }
  values
}

boot_smd_stat <- function(data, indices) {
  d <- data[indices, , drop = FALSE]
  x1 <- d$value[d$cohort == "champagne"]
  x2 <- d$value[d$cohort == "GSE77445"]
  n1 <- length(x1)
  n2 <- length(x2)
  mean_diff <- mean(x1) - mean(x2)
  pooled_sd <- sqrt((((n1 - 1) * var(x1)) + ((n2 - 1) * var(x2))) / (n1 + n2 - 2))
  c(Mean_Difference = mean_diff, Cohens_d = mean_diff / pooled_sd)
}

get_bca_ci <- function(boot_obj, idx, label) {
  ci <- tryCatch(
    boot::boot.ci(boot_obj, index = idx, type = "bca"),
    error = function(e) {
      stop(paste0("BCa CI failed for ", label, ": ", e$message), call. = FALSE)
    }
  )
  if (is.null(ci$bca)) {
    stop(paste0("BCa CI returned no interval for ", label, "."), call. = FALSE)
  }
  ci$bca[4:5]
}

compare_cohorts <- function(champagne_data, champagne_label, gse_data, var_map, R = 10000) {

  result_rows <- lapply(names(var_map), function(champ_col) {

    gse_col <- var_map[[champ_col]]
    label <- paste0(champagne_label, " vs. GSE77445 | ", champ_col)
    message(paste("Bootstrapping cohort difference for:", label))

    x1 <- get_cohort_values(champagne_data, champ_col, champagne_label)
    x2 <- get_cohort_values(gse_data, gse_col, "GSE77445")

    long_df <- data.frame(
      value  = c(x1, x2),
      cohort = factor(
        rep(c("champagne", "GSE77445"), c(length(x1), length(x2))),
        levels = c("champagne", "GSE77445")
      )
    )

    boot_obj <- boot::boot(
      data      = long_df,
      statistic = boot_smd_stat,
      R         = R,
      strata    = long_df$cohort
    )

    ci_mean_diff <- get_bca_ci(boot_obj, 1, paste(label, "(mean difference)"))
    ci_cohens_d  <- get_bca_ci(boot_obj, 2, paste(label, "(Cohen's d)"))

    data.frame(
      Comparison               = paste0(champagne_label, " vs. GSE77445"),
      Variable                 = if (champ_col == "Sex") "Sex (proportion female)" else champ_col,
      Champagne_Column         = champ_col,
      GSE77445_Column          = gse_col,
      N_Champagne              = length(x1),
      N_GSE77445               = length(x2),
      Mean_Champagne           = mean(x1),
      SD_Champagne             = sd(x1),
      Mean_GSE77445            = mean(x2),
      SD_GSE77445              = sd(x2),
      Mean_Difference          = as.numeric(boot_obj$t0[1]),
      Mean_Difference_Boot_SE  = sd(boot_obj$t[, 1]),
      Mean_Difference_Lower_95 = ci_mean_diff[1],
      Mean_Difference_Upper_95 = ci_mean_diff[2],
      Cohens_d                 = as.numeric(boot_obj$t0[2]),
      Cohens_d_Boot_SE         = sd(boot_obj$t[, 2]),
      Cohens_d_Lower_95        = ci_cohens_d[1],
      Cohens_d_Upper_95        = ci_cohens_d[2],
      row.names                = NULL
    )
  })

  do.call(rbind, result_rows)
}

cohort_comparison_results <- bind_rows(
  compare_cohorts(
    champagne_data   = champagne_pre_pred,
    champagne_label  = "champagne_pre",
    gse_data         = GSE77445_pred,
    var_map          = cohort_var_map
  ),
  compare_cohorts(
    champagne_data   = champagne_post_pred,
    champagne_label  = "champagne_post",
    gse_data         = GSE77445_pred,
    var_map          = cohort_var_map
  )
)

cohort_comparison_results$Significant_CI <-
  (cohort_comparison_results$Cohens_d_Lower_95 > 0 & cohort_comparison_results$Cohens_d_Upper_95 > 0) |
  (cohort_comparison_results$Cohens_d_Lower_95 < 0 & cohort_comparison_results$Cohens_d_Upper_95 < 0)

write.csv(
  cohort_comparison_results,
  file.path(main_out_dir, "cohort_comparisons.csv"),
  row.names = FALSE
)