#!/usr/bin/env Rscript
args=(commandArgs(TRUE))

project <- args[1]
target <- args[2]    

library(hipathia)

exp_data <- readRDS(paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/01_mechanistics_models/01_1_mechanistics_models_tcga/",project,"/exp_data_",project,".rds",sep=""))
exp_data[target,] <- exp_data[target,]/0.01	

pathways <- load_pathways("hsa")
results <- hipathia(exp_data, pathways, decompose = FALSE, verbose=FALSE)

path_vals <- get_paths_data(results, matrix = TRUE)
path_normalized <- normalize_paths(path_vals, pathways)

	
file_name <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/04_mm_ko_drug_used/04_3_pathvals_ko_drug_used/",project,"/pathvals_", project,"_",target,".rds",sep="")
	
saveRDS(path_normalized, file=file_name)
