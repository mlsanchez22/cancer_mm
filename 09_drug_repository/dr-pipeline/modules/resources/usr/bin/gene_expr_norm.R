#!/usr/bin/env Rscript

# Purpose: normalize gene expression data with the algorithm of TMM.
# Input 1: path to raw gene counts table, tab-separated format
# Input 2: Boolean, low-expression gene filtering (optional; default FALSE)
# Input 3: Boolean, log2 transformation of gene expression (optional; default TRUE)
# Output 1: TMM normalized gene expression in log scale (rds)
# Output 2: transposed TMM normalized gene expression in log scale (tsv)

library(edgeR)
suppressMessages(library(tibble))


# get arguments from command line

args <- commandArgs(trailingOnly = TRUE)

# check if there are at least one argument, otherwise return an error
if (length(args) < 1) {
  stop("Path to expression data required.", call. = F)
} else if (length(args) == 1) {
  args[2] <- FALSE
  args[3] <- TRUE
} else if (length(args) == 2) {
  args[3] <- TRUE
}

# load resources
message("Loading gene expression...")

expr_file <- args[1]
filt_genes <- as.logical(args[2])
do_log <- as.logical(args[3])

if (!file.exists(expr_file)) {stop("Expression file does not exist.")}
expr <- read.delim(expr_file, row.names = 1, check.names = F)

# read header of input file separately
file_head <- readLines(expr_file, 1)
file_head <- strsplit(file_head, "\t")[[1]]

# direct normalization with no grouping factors
message("Normalizing gene expression...")

expr.TMM <- DGEList(counts = expr)
if (filt_genes) {
  keep <- filterByExpr(expr.TMM)
  expr.TMM <- expr.TMM[keep,,keep.lib.sizes=FALSE]
  message("Keeping ", nrow(expr.TMM), " genes from ", nrow(expr), " after filtering.")
}
expr.TMM <- calcNormFactors(expr.TMM,  method = "TMM")

expr.norm <- cpm(expr.TMM, log = do_log, prior.count = 3)

# write output to tsv
## if colnames have been modified on the process, assign the ones on input file
message("Writing normalized data to table...")

start_colnames <- length(file_head) - length(colnames(expr.norm)) + 1
if (! all.equal(file_head[start_colnames:length(file_head)], colnames(expr.norm))){
  colnames(expr.norm) <- file_head[start_colnames:length(file_head)]
}
saveRDS(expr.norm , file = "tmm-normalized_expr.rds")

expr.norm_tt <- rownames_to_column(as.data.frame(t(expr.norm)), var = "index")
write.table(expr.norm_tt, file = "tmm-normalized_expr_transposed.tsv", sep = "\t", quote = F, row.names=FALSE)


message("Expression data normalization done.")
