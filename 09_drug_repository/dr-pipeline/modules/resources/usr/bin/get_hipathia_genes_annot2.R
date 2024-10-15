#!/usr/bin/env Rscript

# Purpose: add annotations to translated circuit-gene relations.
# Input 1: translated table from previous step of circuit information extraction 
# Input 2, optional: hallmarks table with circuit and hallmark columns
# Output 1: tsv file with extended pathway/circuit/gene(/hallmark) information

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))


# get arguments from command line

args <- commandArgs(trailingOnly = TRUE)

# check if there are at least two argument, otherwise return an error
if (length(args) < 1) {
  stop("Translated hipathia annotations required.", call. = F)
}

# load circuit-genes files

message("Loading pathway resources...")

circ_genes_file <- args[1]
if (!file.exists(circ_genes_file)) stop("Cannot access annotations file.")
circ_genes <- read_tsv(circ_genes_file, show_col_types = F)


### add column for single/indicated circuit effector(s) ###
# Genes indicated in circuit names are used as reference.

message("Establish effector genes...")

single_eff <- circ_genes %>% 
  dplyr::relocate(gene_symbol, .after = "gene") %>% 
  group_by(circuit, node) %>% mutate(id = seq_along(node)) %>% ungroup %>% 
  mutate(eff_node = if_else(eff_node == TRUE, gene_symbol, NA)) %>% 
  tidyr::separate(circuit_name, c(NA,"effector"), sep = ": ", remove = FALSE) %>% 
  tidyr::separate_rows(effector, sep = " ") %>% 
  dplyr::mutate(effector = stringr::str_extract(effector, "[a-zA-Z0-9]+")) %>% 
  group_by(circuit, node) %>% dplyr::mutate(geneiseff = eff_node %in% effector) %>% 
  ungroup %>% 
  dplyr::mutate(eff_gene = if_else(geneiseff == TRUE, eff_node, 
                                   if_else(!is.na(eff_node) & id == 1, eff_node, NA))) %>% 
  dplyr::select(-id, -effector, -geneiseff) %>% distinct()


### if provided, add info on hallmarks/modules ###
# input: file with 2 columns, circuit IDs (circuit) and hallmarks/modules (hallmark)

if (length(args) == 2) {
  message("Adding module/hallmark to circuit information...")
  hallmark_file <- args[2]
  circ_hall <- read_tsv(hallmark_file, show_col_types = F) %>% 
    dplyr::select(circuit, hallmark) %>% distinct %>% 
    summarise(hallmark = paste(hallmark, collapse = ","), .by = "circuit")
  circ_genes <- single_eff %>% 
    left_join(circ_hall, by = "circuit") %>% relocate(hallmark, .after = circuit_name)
} else {
  circ_genes <- single_eff
}

### write the final table of hipathia circuit information ###

message("Writing extended pathway information...")
out_name <- gsub(".tsv$", "-annot.tsv", basename(circ_genes_file))
write_tsv(circ_genes, out_name, na = "")

message("Annotation of hipathia pathways done.")
