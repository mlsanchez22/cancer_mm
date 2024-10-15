#!/usr/bin/env Rscript

# Purpose: creation of a list of relevant genes associated to a given disease.
# Input 1: local path or web link to DISGENET table
# Input 2: disease ID
# Input 3: percentile in range 0-1 above which select highest scoring genes. Default 0.75
# Input 4: order priority by GDA score ("gda") or relevance score ("relevance"; default)
# Output: a list of genes in entrez ID

suppressMessages(library(readr))
suppressMessages(library(dplyr))

# get arguments from command line

args <- commandArgs(trailingOnly = TRUE)

# check if there are at least two arguments, otherwise return an error
if (length(args) < 2) {
  stop("At least a path to DISGENET table and a disease ID must be provided.", call. = F)
} else if (length(args) == 2) {
  # default values
  args[3] <- 0.75
  args[4] <- "relevance"
} else if (length(args) == 3) {
  args[4] <- "relevance"
}

dg_path <- args[1]
disease_id <- args[2]
threshold <- as.numeric(args[3])
order_type <- args[4]

if (findInterval(threshold, c(0,1)) != 1L) {stop("Threshold must be between 0 and 1.")}
if (! order_type %in% c("gda", "relevance")) {stop("Order must be one of gda or relevance.")}

message("Processing data...")

# read table disgenet

dg_table <- read_tsv(dg_path, show_col_types = F)

# check that the disease id is in the table or stop

if (!disease_id %in% dg_table$diseaseId) stop("The disease ID is not valid.", call. = F)

# Filter the table for the requested disease and generate the combined score.

dg_dis <- dplyr::filter(dg_table, diseaseId == disease_id)
if (order_type == "relevance") {
  dg_dis <- dg_dis %>% 
    dplyr::mutate(score = score / DPI * DSI * EI) %>% 
    dplyr::arrange(desc(score))
}

min_sc <- quantile(dg_dis$score, threshold, na.rm = T)
dg_dis.threshold <- filter(dg_dis, score > min_sc)

message(paste("Selected", length(dg_dis.threshold$gene_symbol), 
              "highest scoring genes from", length(dg_dis$geneSymbol), 
              "total genes associated to", dg_dis$diseaseName[1], 
              "by", order_type, "score"))

# write gene list

write_tsv(dplyr::select(dg_dis.threshold, gene_symbol), "geneList.tsv", col_names = F)

message("Gene selection from DisGeNET data done.")
