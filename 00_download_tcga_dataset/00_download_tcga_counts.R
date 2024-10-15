args=(commandArgs(TRUE))

i <- as.numeric(args[1])

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

projects <- getGDCprojects()
projects_tcga <- projects$id[grep("TCGA",projects$id)]

download_tcga_counts <- function(x){

query.exp.hg38 <- GDCquery(project = x,
data.category = "Transcriptome Profiling",
data.type = "Gene Expression Quantification",
workflow.type = "STAR - Counts")

GDCdownload(query.exp.hg38)
expdat <- GDCprepare(
      query = query.exp.hg38,
)

	return(expdat)
}


expdata <- download_tcga_counts(projects_tcga[i])

file <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/00_download_tcga_dataset/",gsub("-","_",tolower(projects_tcga[i])),"/counts_",gsub("-","_",tolower(projects_tcga[i])),".rds",sep="")
saveRDS(expdata, file=file)
