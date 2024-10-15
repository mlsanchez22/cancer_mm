args=(commandArgs(TRUE))

i <- as.numeric(args[1])


library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(hipathia)

projects <- getGDCprojects()
projects_tcga <- projects$id[grep("TCGA",projects$id)]

### patients with drug data

drug <- readRDS(paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_sex/data/",gsub("-","_",tolower(projects_tcga[i])),"/drug_",gsub("-","_",tolower(projects_tcga[i])),".rds",sep=""))
drug$drug_name <- tolower(drug$drug_name)


### patients with clinical data
clinical_patient <- readRDS(paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_sex/data/",gsub("-","_",tolower(projects_tcga[i])),"/clinical_patient_",gsub("-","_",tolower(projects_tcga[i])),".rds",sep=""))
patients_with_clinical_patient_drug <- intersect(clinical_patient$bcr_patient_barcode,drug$bcr_patient_barcode)

sel_clinical_patients <- unique(clinical_patient[unlist(lapply(patients_with_clinical_patient_drug, FUN=function(x) grep(x, clinical_patient$bcr_patient_barcode))),c("bcr_patient_barcode","days_to_birth","gender","vital_status","days_to_death","days_to_last_followup","has_new_tumor_events_information")])
sel_drug <- unique(drug[unlist(lapply(patients_with_clinical_patient_drug, FUN=function(x) grep(x, drug$bcr_patient_barcode))),c("bcr_patient_barcode","drug_name")])

patients_clinical_data_drug <- merge(sel_clinical_patients, sel_drug, by.x="bcr_patient_barcode",by.y="bcr_patient_barcode", all.x=T,all.y=T)
patients_clinical_data_drug$age <- abs(patients_clinical_data_drug$days_to_birth / 365)

### patients with gene expression

counts <- assay(readRDS(paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_sex/data/",gsub("-","_",tolower(projects_tcga[i])),"/counts_",gsub("-","_",tolower(projects_tcga[i])),".rds",sep="")))

id_patients_counts <- unique(
paste("TCGA", 
unlist(lapply(strsplit(colnames(counts),"-"), FUN=function(x) x[2])), 
unlist(lapply(strsplit(colnames(counts),"-"), FUN=function(x) x[3])), sep="-"))

patients_geneexp_clind_drug <- patients_clinical_data_drug[patients_clinical_data_drug$bcr_patient_barcode %in% id_patients_counts,]

saveRDS(patients_geneexp_clind_drug, file=paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/02_selection_samples/02_1_ext_information_patients_tcga/",
gsub("-","_",tolower(projects_tcga[i])),"/patients_geneexp_clind_drug_",gsub("-","_",tolower(projects_tcga[i])),".rds",sep=""))


### sample_type

sample_type <- unlist(lapply(strsplit(colnames(counts),"-"), FUN=function(x) x[4]))

sample_type[grep("11", sample_type)] <- "Solid Tissue Normal"
sample_type[grep("01", sample_type)] <- "Primary Solid Tumor"
sample_type[grep("02", sample_type)] <- "Recurrent Solid Tumor"
sample_type[grep("03", sample_type)] <- "Primary Blood Derived Cancer"
sample_type[grep("04", sample_type)] <- "Recurrent Blood Derived Cancer"
sample_type[grep("05", sample_type)] <- "Primary Solid Tumor"
sample_type[grep("06", sample_type)] <- "Metastatic"
sample_type[grep("07", sample_type)] <- "Metastatic"

df_patient_sample_type <- data.frame("patient" = colnames(counts), "sample_type"=sample_type)

df_patient_sample_type$pat <- paste(unlist(lapply(df_patient_sample_type$patient, FUN=function(x) unlist(strsplit(x,"-"))[1])), unlist(lapply(df_patient_sample_type$patient, FUN=function(x) unlist(strsplit(x,"-"))[2])),unlist(lapply(df_patient_sample_type$patient, FUN=function(x) unlist(strsplit(x,"-"))[3])), sep="-")

patients_no_use <- unlist(lapply(as.data.frame(df_patient_sample_type %>% filter(sample_type == "Primary Solid Tumor") %>% group_by(pat) %>% summarise(n=n())  %>% filter(n!=1) %>% select(pat))$pat, FUN=function(x) names(apply(counts[,df_patient_sample_type[df_patient_sample_type$pat %in% x,"patient"]], 2, var))[apply(counts[,df_patient_sample_type[df_patient_sample_type$pat %in% x,"patient"]], 2, var) != min(apply(counts[,df_patient_sample_type[df_patient_sample_type$pat %in% x,"patient"]], 2, var))]))

patient_sample_type <- df_patient_sample_type[!(df_patient_sample_type$patient %in% patients_no_use),]

patients_geneexp_clind_drug_samtype <- merge(patients_geneexp_clind_drug, patient_sample_type, by.x="bcr_patient_barcode", by.y="pat")

saveRDS(patients_geneexp_clind_drug_samtype, file=paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/02_selection_samples/02_1_ext_information_patients_tcga/",
gsub("-","_",tolower(projects_tcga[i])),"/patients_geneexp_clind_drug_samptype_",gsub("-","_",tolower(projects_tcga[i])),".rds",sep=""))


