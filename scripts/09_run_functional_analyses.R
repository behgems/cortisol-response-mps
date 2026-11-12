# 09_run_functional_analyses.R
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(MethylToSNP)
library(openxlsx)
library(missMethyl)
library(future)
library(furrr)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(readr)
library(forcats)
library(ggtext)

ens <- readRDS("./output/02_generate_mps/ensemble/ens_cortisol_auc.rds")
gc_sensitive_sites <- openxlsx::read.xlsx("./data/provencal_2020/pnas_1820842116_sd01_S7.xlsx")
epistress_sites <- openxlsx::read.xlsx("./data/provencal_2020/pnas_1820842116_sd01_S8.xlsx")
mcig_probes_path <- "./data/mcigarette/Predictors_long.csv"
required_probes_df <- read.csv("./output/04_analyse_shap_importance/shap_global_importance_comparison.csv")
champagne_pre_proc <- readRDS("./output/03_generate_predictions/champagne_pre_pred/ens_cortisol_auc_on_champagne_pre_pred_proc.rds")
champagne_post_proc <- readRDS("./output/03_generate_predictions/champagne_post_pred/ens_cortisol_auc_on_champagne_post_pred_proc.rds")
GSE77445_proc <- readRDS("./output/02_generate_mps/ensemble/metamodel_processing/01_process_data/processed_data.rds")
GSE77445_raw <- readRDS("./output/01_prepare_data/GSE77445_formatted_beta.rds")

genes_of_interest <- c(
  "CRH", "CRHR1", "FKBP5", "HSD11B1", "HSD11B2", "LEP", "MC2R", "NR3C1", "NR3C2", "OXTR", "SKA2", "SLC6A4"
)

n_cores <- 10
stress_responsiveness_alpha <- 0.05
main_out_dir <- file.path("./output/09_run_functional_analyses")

if (!dir.exists(main_out_dir)) dir.create(main_out_dir, recursive = TRUE)

set.seed(315)

plan(multisession, workers = n_cores)

# Get annotation
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

# Add annotation to importance df
names(required_probes_df)[names(required_probes_df) == "probe"] <- "Name"
required_probes_df <- merge(required_probes_df, anno, by = "Name")
write.csv(required_probes_df, file.path(main_out_dir, "required_probes_anno.csv"), row.names = FALSE)

# Get required probes vector
required_probes_vec <- required_probes_df$Name

# Get probes that had positive importance in any of the data sets
positive_probes <- required_probes_df[required_probes_df$mean_scaled_importance > 0, ]
write.csv(positive_probes, file.path(main_out_dir, "positive_probes_anno.csv"), row.names = FALSE)
positive_probes <- positive_probes$Name

# Get background probes from full training data set and add anno
background_probes <- rownames(GSE77445_proc$methyl)
background_probes <- data.frame(Name = background_probes)
background_probes_anno <- anno[anno$Name %in% background_probes$Name,]                                        
background_probes <- merge(background_probes, background_probes_anno, by = "Name")
write.csv(background_probes, file.path(main_out_dir, "background_probes_anno.csv"), row.names = FALSE)
background_probes <- background_probes$Name

message("Checking for probes in the training set that might be biased by SNPs, based on their patterns of methylation...")
snp_info <- MethylToSNP::MethylToSNP(GSE77445_raw$methyl, verbose = TRUE)
print(snp_info)    

message("Checking the number of probes that overlap with those previously identified as GC-sensitive...")
sink(file.path(main_out_dir, "gc_probe_comparison.txt"), append = FALSE)

print("Total GC-sensitive and EpiStress sites:")
print(length(unique(c(gc_sensitive_sites$CpG, epistress_sites$CpG))))

gc_sensitive_sites_overlap <- required_probes_vec[required_probes_vec %in% gc_sensitive_sites$CpG]
epistress_sites_overlap <- required_probes_vec[required_probes_vec %in% epistress_sites$CpG]

print("Number of required probes that were previously identified as GC-sensitive:")
print(length(gc_sensitive_sites_overlap))
print("Probe IDs:")
print(gc_sensitive_sites_overlap)

print("Number of required probes that are in EpiStress:")
print(length(epistress_sites_overlap))
print("Probe IDs:")
print(epistress_sites_overlap)
sink()

# Relation to island permutation tests
permutation_test_for_relation_to_island <- function(
  category_to_test,
  n_permutations,
  background_probes,
  model_probes,
  cg_anno_df
) {
  probes_in_category <- cg_anno_df$Name[cg_anno_df$Relation_to_Island == category_to_test]
  
  observed_proportion <- sum(model_probes %in% probes_in_category) / length(model_probes)
  null_proportion <- sum(background_probes %in% probes_in_category) / length(background_probes)
  
  permutation_proportions <- furrr::future_map_dbl(1:n_permutations, function(i) {
    random_sample <- sample(background_probes, size = length(model_probes))
    sum(random_sample %in% probes_in_category) / length(random_sample)
  }, .options = furrr_options(seed = TRUE))
  
  p_value_overrepresentation <- (sum(permutation_proportions >= observed_proportion) + 1) / (n_permutations + 1)
  p_value_underrepresentation <- (sum(permutation_proportions <= observed_proportion) + 1) / (n_permutations + 1)
    
  mean_permuted <- mean(permutation_proportions)
  abs_diff_observed <- abs(observed_proportion - mean_permuted)
  abs_diff_permuted <- abs(permutation_proportions - mean_permuted)
  p_value_two_sided <- (sum(abs_diff_permuted >= abs_diff_observed) + 1) / (n_permutations + 1)
  
  list(
    observed_proportion = observed_proportion,
    null_proportion = null_proportion,
    mean_permuted_proportion = mean_permuted,
    p_value_overrepresentation = p_value_overrepresentation,
    p_value_underrepresentation = p_value_underrepresentation,
    p_value_two_sided = p_value_two_sided
  )
}
                          
# Check model probes for genomic feature annotations
permutation_test_for_genomic_feature_annotations <- function(
  category_to_test,
  n_permutations,
  background_probes,
  model_probes,
  cg_anno_df
) {
    is_category_present <- sapply(strsplit(cg_anno_df$UCSC_RefGene_Group, ";"), function(x) category_to_test %in% x)
    probes_in_category <- cg_anno_df$Name[is_category_present]
          
    observed_proportion <- sum(model_probes %in% probes_in_category) / length(model_probes)
    null_proportion <- sum(background_probes %in% probes_in_category) / length(background_probes)
  
    permutation_proportions <- furrr::future_map_dbl(1:n_permutations, function(i) {
      random_sample <- sample(background_probes, size = length(model_probes))
      sum(random_sample %in% probes_in_category) / length(random_sample)
    }, .options = furrr_options(seed = TRUE))
  
  p_value_overrepresentation <- (sum(permutation_proportions >= observed_proportion) + 1) / (n_permutations + 1)
  p_value_underrepresentation <- (sum(permutation_proportions <= observed_proportion) + 1) / (n_permutations + 1)
  
  mean_permuted <- mean(permutation_proportions)
  abs_diff_observed <- abs(observed_proportion - mean_permuted)
  abs_diff_permuted <- abs(permutation_proportions - mean_permuted)
  p_value_two_sided <- (sum(abs_diff_permuted >= abs_diff_observed) + 1) / (n_permutations + 1)
  
  list(
    observed_proportion = observed_proportion,
    null_proportion = null_proportion,
    mean_permuted_proportion = mean_permuted,
    p_value_overrepresentation = p_value_overrepresentation,
    p_value_underrepresentation = p_value_underrepresentation,
    p_value_two_sided = p_value_two_sided
  )
}

# Perform permutation tests for positive (i.e., important) probes
message("Performing relation to island and genomic feature permutation tests for important probes...")
sink(file.path(main_out_dir, "permutation_tests_important_probes.txt"), append = FALSE)

# Relation to island tests
for (category in c("Island", "N_Shore", "S_Shore", "OpenSea")) {
    result <- permutation_test_for_relation_to_island(
        category_to_test = category,
        n_permutations = 1000,
        background_probes = background_probes,
        model_probes = positive_probes,
        cg_anno_df = anno
    )
    print(paste("Category:", category))
    print(result)
}

# Genomic feature tests
for (category in c("3'UTR", "5'UTR", "1stExon", "Body", "TSS1500", "TSS200")) {
    result <- permutation_test_for_genomic_feature_annotations(
        category_to_test = category,
        n_permutations = 1000,
        background_probes = background_probes,
        model_probes = positive_probes,
        cg_anno_df = anno
    )
    print(paste("Category:", category))
    print(result)
}
sink()

message("Performing GO and KEGG ORA for important probes...")
tasks <- list(
  list(
    sig_probes = positive_probes,
    collection = "GO",
    output_filename = "go_results_important_probes.csv"
  ),
  list(
    sig_probes = positive_probes,
    collection = "KEGG",
    output_filename = "kegg_results_important_probes.csv"
  )
)

run_gsea_task <- function(task_params) {  
  message(paste("Starting analysis for:", task_params$output_filename))
  gsea_results <- missMethyl::gometh(
    sig.cpg = task_params$sig_probes,
    all.cpg = background_probes,
    collection = task_params$collection,
    array.type = "450K",
    anno = anno,
    sig.genes = TRUE
  )
  gsea_results <- gsea_results[order(gsea_results$P.DE, decreasing = FALSE), ]  
  output_path <- file.path(main_out_dir, task_params$output_filename)
  write.csv(gsea_results, output_path, row.names = FALSE)
  message(paste("Finished and saved:", task_params$output_filename))
  return(output_path)
}

furrr::future_walk(
  .x = tasks, 
  .f = ~run_gsea_task(.x),
  .options = furrr_options(
      globals = c("run_gsea_task", "background_probes", "anno", "main_out_dir"), 
      seed = TRUE
  )
)

plan(sequential)

message("Checking if the required probes are annotated to genes of a priori interest...")
valid_gene_names <- !is.na(required_probes_df$UCSC_RefGene_Name) & required_probes_df$UCSC_RefGene_Name != ""
all_genes_list <- strsplit(required_probes_df$UCSC_RefGene_Name[valid_gene_names], ";")
unique_genes <- unique(unlist(all_genes_list))
hits <- genes_of_interest[genes_of_interest %in% unique_genes]

if (length(hits) > 0) {
    search_pattern <- paste(hits, collapse = "|")
    matching_rows <- grepl(search_pattern, required_probes_df$UCSC_RefGene_Name)
    filtered_df <- required_probes_df[matching_rows, ]
    write.csv(filtered_df, file.path(main_out_dir, "genes_of_interest.csv"), row.names = FALSE)
}

message("Performing paired t-tests for pre- vs. post-stress methylation change in Champagne data...")

# Check sample alignment
if (!all(rownames(champagne_pre_proc$pheno) == colnames(champagne_pre_proc$methyl))) {
  stop("Pre-stress pheno and methyl data are not aligned")
}
if (!all(rownames(champagne_post_proc$pheno) == colnames(champagne_post_proc$methyl))) {
  stop("Post-stress pheno and methyl data are not aligned")
}

# Set methylation colnames to participant IDs rather than sample IDs
colnames(champagne_pre_proc$methyl) <- champagne_pre_proc$pheno$PID
colnames(champagne_post_proc$methyl) <- champagne_post_proc$pheno$PID

# Get the methylation matrices
pre_methyl <- champagne_pre_proc$methyl
post_methyl <- champagne_post_proc$methyl

# Confirm that both matrices have the same probes in the same order
common_probes <- intersect(rownames(pre_methyl), rownames(post_methyl))
pre_methyl <- pre_methyl[common_probes, ]
post_methyl <- post_methyl[common_probes, ]

# Ensure both matrices have the same samples in the same order
common_samples <- intersect(colnames(pre_methyl), colnames(post_methyl))
pre_methyl <- pre_methyl[, common_samples]
post_methyl <- post_methyl[, common_samples]

# Check if data are valid
if (!all(rownames(pre_methyl) == rownames(post_methyl)) || !all(colnames(pre_methyl) == colnames(post_methyl))) {
  stop("Error: pre-stress and post-stress methylation matrices do not have matching dimensions or order after filtering")
}
if (length(common_samples) < 2) {
    stop("Error: fewer than two common samples between pre-stress and post-stress data sets")
}

ttest_results_file <- file.path(main_out_dir, "pre_post_stress_ttest_results.csv")

message("Running t-tests...")

# Run t-tests
fast_t_tests <- matrixTests::row_t_paired(x = post_methyl, y = pre_methyl)

ttest_results <- data.frame(
probe_id = rownames(fast_t_tests),
t_statistic = fast_t_tests$statistic,
p_value = fast_t_tests$pvalue,
mean_diff = fast_t_tests$mean.diff,
stringsAsFactors = FALSE
)

# Remove row names and sort by p value
rownames(ttest_results) <- NULL
ttest_results <- ttest_results[order(ttest_results$p_value), ]

write.csv(
ttest_results,
ttest_results_file,
row.names = FALSE
)
message(paste("T-tests complete. Results saved to:", ttest_results_file))

message("Performing enrichment analysis for stress-responsive probes in Champagne data...")
gse77445_background <- rownames(GSE77445_proc$methyl) # All probes available for model training in GSE77445
champagne_background <- rownames(pre_methyl) # All probes available for the t-tests in the Champagne data
final_background_probes <- intersect(gse77445_background, champagne_background) # Background is the intersection of the two

# Filter the probe sets to the background
positive_probes_in_common <- intersect(positive_probes, final_background_probes)
required_probes_in_common <- intersect(required_probes_vec, final_background_probes)

# Identify nominally stress-responsive probes from the t-test results
significant_probes <- ttest_results$probe_id[which(ttest_results$p_value < stress_responsiveness_alpha)]

permutation_enrichment_test <- function(probes_to_test, significant_probes, background_probes, n_permutations = 1000) {
    observed_overlap <- sum(probes_to_test %in% significant_probes)
    n_to_sample <- length(probes_to_test)
    
    cat("Running", n_permutations, "permutations for:", deparse(substitute(probes_to_test)))
    
    # Run in parallel
    permuted_overlaps <- furrr::future_map_dbl(1:n_permutations, function(i) {
        random_sample <- sample(background_probes, size = n_to_sample, replace = FALSE)
        sum(random_sample %in% significant_probes)
    }, .options = furrr_options(seed = TRUE))
    
    p_value_enrichment <- (sum(permuted_overlaps >= observed_overlap) + 1) / (n_permutations + 1)
    p_value_depletion <- (sum(permuted_overlaps <= observed_overlap) + 1) / (n_permutations + 1)
    
    mean_permuted <- mean(permuted_overlaps)
    abs_diff_observed <- abs(observed_overlap - mean_permuted)
    abs_diff_permuted <- abs(permuted_overlaps - mean_permuted)
    p_value_two_sided <- (sum(abs_diff_permuted >= abs_diff_observed) + 1) / (n_permutations + 1)
    
    list(
        probeset_name = deparse(substitute(probes_to_test)),
        observed_overlap = observed_overlap,
        mean_permuted_overlap = mean_permuted,
        n_permutations = n_permutations,
        p_value_enrichment = p_value_enrichment,
        p_value_depletion = p_value_depletion,
        p_value_two_sided = p_value_two_sided
    )
}

# Run stress responsive probe ORA
sink(file.path(main_out_dir, "stress_responsive_enrichment_important_probes.txt"))
perm_results_positive <- permutation_enrichment_test(
    probes_to_test = positive_probes_in_common,
    significant_probes = significant_probes,
    background_probes = final_background_probes,
    n_permutations = 1000
)
cat("\nPermutation test results:\n")
cat("Observed overlap:", perm_results_positive$observed_overlap, "\n")
cat("Mean permuted overlap:", round(perm_results_positive$mean_permuted_overlap, 2), "\n")
cat("p value for enrichment (one-sided):", perm_results_positive$p_value_enrichment, "\n")
cat("p value for depletion (one-sided):", perm_results_positive$p_value_depletion, "\n")
cat("p value (two-sided):", perm_results_positive$p_value_two_sided, "\n")
sink()