# 01_prepare_data.R
library(GEOquery)
library(Biobase)

champagne_beta_path <- file.path("./data/methyl/champagne/beta_no_gap_probes.rda")
champagne_pheno_path <- file.path("./data/methyl/champagne/DataAcuteStressMeth.xlsx")
GSE77445_path <- file.path("./data/methyl/GSE77445/GSE77445_series_matrix.txt")
GSE119078_path <- file.path("./data/methyl/GSE119078/GSE119078_series_matrix.txt.gz")
main_out_dir <- file.path("./output/01_prepare_data")
dir.create(main_out_dir, showWarnings = FALSE, recursive = TRUE)

message("Preparing GSE119078 methylation data...")
gse <- GEOquery::getGEO(filename = GSE119078_path, getGPL = FALSE)
pheno <- Biobase::pData(gse) 
methyl <- Biobase::exprs(gse)
common_samples <- intersect(colnames(methyl), rownames(pheno))
pheno <- pheno[common_samples, , drop = FALSE]
methyl <- methyl[, common_samples]
pheno <- pheno[pheno[["disease status:ch1"]] == "Control", , drop = FALSE]
methyl <- methyl[, rownames(pheno)]

GSE119078_formatted_beta <- list(pheno = pheno, methyl = methyl)
saveRDS(GSE119078_formatted_beta, file.path(main_out_dir, "GSE119078_formatted_beta.rds"))

message("Preparing Champagne methylation data...")
load(file = champagne_beta_path)
champagne_pheno <- openxlsx::read.xlsx(champagne_pheno_path)

common_samples <- intersect(colnames(beta_no_gap_probes), champagne_pheno$Basename)
pheno_ordered <- champagne_pheno[match(common_samples, champagne_pheno$Basename), ]
beta_ordered <- beta_no_gap_probes[, common_samples]

combined_data <- list(pheno = as.data.frame(pheno_ordered), methyl = as.matrix(beta_ordered))
rownames(combined_data$pheno) <- colnames(combined_data$methyl)

pheno_names <- c("PID", "Condition", "Age", "Time", "Sex", "Ethnicity_Race", "Basename", "AUCi_scc", "s1")
combined_data$pheno <- combined_data$pheno[, pheno_names, drop = FALSE]

combined_data$pheno$Ethnicity_Race <- as.factor(combined_data$pheno$Ethnicity_Race)
combined_data$pheno$Ethnicity_Race <- relevel(combined_data$pheno$Ethnicity_Race, "White")
combined_data$pheno$Sex <- as.factor(combined_data$pheno$Sex)
combined_data$pheno$Sex <- relevel(combined_data$pheno$Sex, "female")
                              
combined_data$pheno <- combined_data$pheno[combined_data$pheno$Condition == 'E', ]
combined_data$methyl <- combined_data$methyl[, rownames(combined_data$pheno)]

champagne_pre_formatted_beta <- combined_data
champagne_post_formatted_beta <- combined_data

champagne_pre_formatted_beta$pheno <- champagne_pre_formatted_beta$pheno[champagne_pre_formatted_beta$pheno$Time == 'pre', ]
champagne_pre_formatted_beta$methyl <- champagne_pre_formatted_beta$methyl[, rownames(champagne_pre_formatted_beta$pheno)]

champagne_post_formatted_beta$pheno <- champagne_post_formatted_beta$pheno[champagne_post_formatted_beta$pheno$Time == 'post', ]
champagne_post_formatted_beta$methyl <- champagne_post_formatted_beta$methyl[, rownames(champagne_post_formatted_beta$pheno)]

saveRDS(champagne_pre_formatted_beta, file.path(main_out_dir, "champagne_pre_formatted_beta.rds"))
saveRDS(champagne_post_formatted_beta, file.path(main_out_dir, "champagne_post_formatted_beta.rds"))

message("Preparing GSE77445 methylation data...")
gse <- GEOquery::getGEO(filename = GSE77445_path, getGPL = FALSE)
pheno_data <- Biobase::pData(gse) 
pheno <- pheno_data[, c(
    "cortisol stress response area under the curve (auc) with respect to the increase:ch1",
    "age:ch1", "gender:ch1", "total score on the childhood trauma questionnaire:ch1",
    "total score on the life stressor checklist- revised:ch1"
), drop = FALSE]

colnames(pheno) <- c("cortisol_auc", "age", "gender", "childhood_trauma", "life_stress")

for (col_name in names(pheno)) {
    pheno[[col_name]][pheno[[col_name]] == "--"] <- NA
    pheno[[col_name]][pheno[[col_name]] == "NA"] <- NA
    if (col_name %in% c("gender")) {
        pheno[[col_name]] <- as.factor(pheno[[col_name]])
    } else {
        pheno[[col_name]] <- as.numeric(pheno[[col_name]])
    }
}

pheno$gender <- relevel(pheno$gender, "Female")
methyl <- Biobase::exprs(gse)
common_samples <- intersect(colnames(methyl), rownames(pheno))
pheno <- pheno[common_samples, , drop = FALSE]
methyl <- methyl[, common_samples]

GSE77445_formatted_beta <- list(pheno = pheno, methyl = methyl)
saveRDS(GSE77445_formatted_beta, file.path(main_out_dir, "GSE77445_formatted_beta.rds"))