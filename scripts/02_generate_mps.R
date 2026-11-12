# 02_generate_mps.R
library(lgr)
library(caret)
library(dplyr)
library(ranger)
library(infotheo)
library(glmnet)
library(Metrics)
library(doParallel)
library(foreach)
library(doRNG)
library(ggplot2)
library(stringr)

source("./scripts/generate_mps/utils.R")
source("./scripts/generate_mps/00_outer_loop.R")
source("./scripts/generate_mps/01_process_data.R")
source("./scripts/generate_mps/02_rank_probes.R")
source("./scripts/generate_mps/03_select_probes.R")
source("./scripts/generate_mps/04_train_models.R")
source("./scripts/generate_mps/05_test_models.R")
source("./scripts/generate_mps/06_aggregate_results.R")
source("./scripts/generate_mps/07_ensemble_models.R")

raw_data <- readRDS("./output/01_prepare_data/GSE77445_formatted_beta.rds")
config <- readRDS("./data/config.rds")
main_out_dir <- file.path("output/02_generate_mps")

outer_loop(
  config = config,
  data = raw_data,
  out_dir = main_out_dir
)