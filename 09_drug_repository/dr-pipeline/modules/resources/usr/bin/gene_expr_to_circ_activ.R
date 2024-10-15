#!/usr/bin/env Rscript

# Purpose: calculate circuit activity from log-normalized counts.
# Input 1: path to normalized expression data rds
# Input 2: path to hipathia pathways rds
# Input 3: pathway id 
# Output 1: transposed circuit activity values (at pathway-level)

suppressMessages(library(hipathia))
suppressMessages(library(tibble))

# get arguments from command line

args <- commandArgs(trailingOnly = TRUE)

#check if there are at least two arguments, otherwise return an error
if (length(args) < 2) {
  stop("Paths to expression data and pathways rds required.", call. = F)
} 

# load resources

message("Loading normalized expression and pathways resources...")

expression_file <- args[1]
pathways_rds <- args[2]


if (!file.exists(expression_file)) {stop("Normalized expression file does not exist.")}
if (!file.exists(pathways_rds)) {stop("Pathways rds file does not exist.")}
logCPM <- readRDS(expression_file)
pathways <- readRDS(pathways_rds)

# read header of input file separately
file_head <- readLines(expression_file, 1)
file_head <- strsplit(file_head, "\t")[[1]]

# filtering by pathway id
message("Selecting pathway...")
path_id <- args[3]
if (!any(names(pathways$pathigraphs)%in%path_id)) stop("Pathway not found.") # error if the pathway is not contained in the general pathways object
pathways$pathigraphs <- pathways$pathigraphs[path_id]

# run hipathia
message("Computing activation signal...")

logCPM <- as.matrix(logCPM)
results <- hipathia(logCPM, pathways, verbose = TRUE)
path_vals <- get_paths_data(results, matrix = T)
path_vals <- normalize_paths(path_vals, pathways)


# write output
## if colnames have been modified, assign the ones on input file
message("Writing results to RDS file...")

path_vals_tt <- rownames_to_column(as.data.frame(t(path_vals)), var = "index")

saveRDS(path_vals_tt , file = paste0("path_vals_transposed-" , path_id , ".rds"))

message("Computation of subpathway activation level done.")
