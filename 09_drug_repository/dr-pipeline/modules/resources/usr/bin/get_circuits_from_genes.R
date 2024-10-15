#!/usr/bin/env Rscript

# Purpose: from a list of genes, extract the associated circuits.
# Input 1: annotated table of circuit-gene associations in hipathia
# Input 2: a list of genes in either entrez or symbol ID, 1 column, no header
# Input 3: gene column name for the same IDs in the circuit-gene table
# Input 4: circuit column name in the circuit-gene table
# Input 5, optional: minimum number of genes per circuit; defaults to 1
# Output 1: a list of circuits

suppressMessages(library(readr))
suppressMessages(library(dplyr))

# get arguments from command line

args <- commandArgs(trailingOnly = TRUE)

# check if there are at least two arguments, otherwise return an error

if (length(args) < 4) {
  stop("Gene list, hipathia annotations and gene/circuit column names required.", call. = F)
}

# load data
message("Loading data...")

circ_genes_file <- args[1]
if (!file.exists(circ_genes_file)) stop("Cannot access hipathia annotations file.")
circ_genes <- read_tsv(circ_genes_file, show_col_types = F)
head(circ_genes)

list_genes_file <- args[2]
if (!file.exists(list_genes_file)) stop("Cannot access gene list.")
list_genes <- read_tsv(list_genes_file, col_names = "id", show_col_types = F)

if (length(args) == 5) {
  min_gene <- as.numeric(args[5])
} else {
  min_gene <- 1
}

# filter annotations based on gene list and minimum associated genes
message("Selecting circuits from genes...")

filt_annot <- circ_genes %>% 
  dplyr::rename(id = args[3], !!args[4] := args[4]) %>% 
  dplyr::filter(id %in% list_genes$id)

if (nrow(filt_annot) == 0){stop("No genes in circuits.")}

filt_annot <- filt_annot %>% group_by(!!args[4]) %>% 
  dplyr::mutate(n_genes = length(unique(id))) %>% ungroup %>% 
  dplyr::filter(n_genes >= min_gene)

if (nrow(filt_annot) == 0){stop(paste("No circuits with at least", min_gene, "genes per circuit."))}

message(paste("From a list of", length(unique(list_genes$id)), 
              "genes,", length(unique(filt_annot$id)), 
              "are associated to", length(unique(filt_annot[[args[4]]])), 
              "circuits with a minimum of", min_gene, "gene(s) per circuit."))

# export circuits to a table without header

write_tsv(unique(dplyr::select(filt_annot, !!args[4])),
          "circuitList.tsv", col_names = F)

message("Conversion from gene list to circuit list done.")
