#!/usr/bin/env Rscript
args=(commandArgs(TRUE))
x <- args[1]


library(hipathia)
	
exp_data <- readRDS("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/01_mechanistics_models/01_2_mechanistics_models_gtex/exp_data_gtex.rds")
exp_data[x,] <- exp_data[x,]*0.01	

pathways <- load_pathways("hsa")
results <- hipathia(exp_data, pathways, decompose = FALSE, verbose=FALSE)

path_vals <- get_paths_data(results, matrix = TRUE)

path_normalized <- normalize_paths(path_vals, pathways)
	
file_name <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/04_mm_ko_drug_used/04_4_pathvals_ko_drug_used_gtex/pathvals_", x, "_gtex.rds",sep="")
	
saveRDS(path_normalized, file=file_name)
