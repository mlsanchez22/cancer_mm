args=(commandArgs(TRUE))

i <- as.numeric(args[1])

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

projects <- getGDCprojects()
projects_tcga <- projects$id[grep("TCGA",projects$id)]

download_tcga_drug <- function(x){
    
    query <- GDCquery(
        project = x,
        data.category = "Clinical",
        data.type = "Clinical Supplement",
        data.format = "BCR XML"
    )
    GDCdownload(query)
    drug <- GDCprepare_clinic(query,"drug")
    return(drug)
}



drug <- download_tcga_drug(projects_tcga[i])

file <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/00_download_tcga_dataset/",gsub("-","_",tolower(projects_tcga[i])),"/drug_",gsub("-","_",tolower(projects_tcga[i])),".rds",sep="")
saveRDS(drug, file=file)

