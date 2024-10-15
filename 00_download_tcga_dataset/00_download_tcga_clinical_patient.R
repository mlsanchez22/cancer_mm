args=(commandArgs(TRUE))

i <- as.numeric(args[1])

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)

projects <- getGDCprojects()
projects_tcga <- projects$id[grep("TCGA",projects$id)]

download_tcga_clinical_patient <- function(x){

        query <- GDCquery(project = x,
                data.category = "Clinical", 
                file.type = "xml")
                GDCdownload(query)

        clinical_patient <- GDCprepare_clinic(query, clinical.info = "patient")
  
        return(clinical_patient)

}



clinical_patient <- download_tcga_clinical_patient(projects_tcga[i])

file <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/00_download_tcga_dataset/",gsub("-","_",tolower(projects_tcga[i])),"/clinical_patient_",gsub("-","_",tolower(projects_tcga[i])),".rds",sep="")
saveRDS(clinical_patient, file=file)

