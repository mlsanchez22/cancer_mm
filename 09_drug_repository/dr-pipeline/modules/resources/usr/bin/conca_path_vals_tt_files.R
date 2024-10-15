#!/usr/bin/env Rscript

# Purpose: concatenate circuit activity from path vals (tt) rds files.
# Input 1: path to path vals transposed data
# Output 1: circuit activity values

message("Reading circuit activation matrices at pathway-level...")

path_vals_tt_files <- list.files()[grep("path_vals_transposed-" , list.files())]

path_vals_tt <- lapply(path_vals_tt_files, readRDS)

path_vals_tt_df <- Reduce(merge , path_vals_tt)

message("Writing results...")

write.table(path_vals_tt_df, file = "path_vals_transposed_concatenated.tsv", sep = "\t", quote = F, row.names=FALSE)
