#!/usr/bin/env Rscript
args=(commandArgs(TRUE))

project <- args[1]

library(survival)
library(hipathia)
library(dplyr)

pathways <- load_pathways("hsa")

pats_geneexp_clind_drug_samptype_file <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/02_selection_samples/02_1_ext_information_patients_tcga/",project,"/patients_geneexp_clind_drug_samptype_",project,".rds", sep="")

pats_geneexp_clind_drug_samptype <- readRDS(pats_geneexp_clind_drug_samptype_file)

pats_geneexp_clind_drug_samptype$time <- pats_geneexp_clind_drug_samptype$days_to_last_followup
pats_geneexp_clind_drug_samptype$time[is.na(pats_geneexp_clind_drug_samptype$days_to_last_followup)] <- pats_geneexp_clind_drug_samptype$days_to_death[is.na(pats_geneexp_clind_drug_samptype$days_to_last_followup)]

patient_sex_rec_time <- unique(pats_geneexp_clind_drug_samptype %>% filter(sample_type == "Primary Solid Tumor") %>% select(gender, has_new_tumor_events_information, patient,time))

patient_sex_rec_time <- patient_sex_rec_time %>% mutate(bi_has_new_tumor_events_information = case_when(
    has_new_tumor_events_information == "YES" ~ 1,
    has_new_tumor_events_information == "NO" ~ 0))
    
# get pathvals of patients tumor

pathvals_file <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/01_mechanistics_models/01_1_mechanistics_models_tcga/",project,"/pathvals_",project,".rds",sep="")    
pathvals <- readRDS(pathvals_file)

physiological_paths <- read.table("/mnt/lustre/scratch/CBRA/projects/heterogeneity/data/physiological_paths.tsv", sep="\t") 

pathvals_tumor <- pathvals[unlist(lapply(physiological_paths$V2, FUN=function(x) grep(x, rownames(pathvals)))),as.data.frame(pats_geneexp_clind_drug_samptype %>% filter(sample_type == "Primary Solid Tumor"))$patient]

# scale values of pathvals selecting +-0.5

pathvals_tumor_sel_scale <- t(pathvals_tumor[rownames(pathvals_tumor)[apply(abs(scale(pathvals_tumor)) > 0.5, 1, sum) == length(colnames(pathvals_tumor))],])
colnames(pathvals_tumor_sel_scale) <- get_path_names(pathways, colnames(pathvals_tumor_sel_scale))

pat_id <- rownames(pathvals_tumor_sel_scale)
pathvals_tumor_sel_scale <- as.data.frame(pathvals_tumor_sel_scale)
pathvals_tumor_sel_scale$patients <- pat_id


### make df for the calculations

patient_sex_rec_time_circ_scale <- merge(patient_sex_rec_time,pathvals_tumor_sel_scale, by.x="patient", by.y="patients")
pathvals_tumor_sel_scale$patients <- NULL

circs_fscale <- colnames(pathvals_tumor_sel_scale)

### Cox PH model

# Function to compute the p-value of the interaction term in Cox PH model with SEX 

pvalue_coxph_circ_sex <- function(x) {
	
	# Construct the formula as a string and then convert it to a formula object
	formula_str <- paste0("Surv(time, bi_has_new_tumor_events_information) ~ gender*`", x, "`")
	coxph_model <- coxph(as.formula(formula_str), data = patient_sex_rec_time_circ_scale)
	cox_summary <- summary(coxph_model)
	pvalue <- cox_summary$coefficients[3,5]
	return(pvalue)

}                  

pvalues_sex_coxph <- lapply(circs_fscale, FUN = function(x) pvalue_coxph_circ_sex(x))

# Name the p-values with the respective column names for easier interpretation

names(pvalues_sex_coxph) <- circs_fscale

## adjust pv

pvadj_sex_coxph <- p.adjust(as.numeric(unlist(pvalues_sex_coxph)), method="fdr")      
names(pvadj_sex_coxph) <- circs_fscale                

save_file_coxph_sex <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/07_survival_analysis/",project,"/coxph_pvadj_sex_",project,".rds",sep="")

saveRDS(pvadj_sex_coxph, file=save_file_coxph_sex)              

circ_sex_coxph <- na.omit(names(pvadj_sex_coxph)[pvadj_sex_coxph < 0.05])


if (length(circ_sex_coxph) > 0){
    
pvalue_zph_test_circ_sex <- function(x) {

formula_str <- paste0("Surv(time, bi_has_new_tumor_events_information) ~ gender * `", x, "`")
  coxph_model <- coxph(as.formula(formula_str), data = patient_sex_rec_time_circ_scale)
  tryCatch({
    zph_test <- cox.zph(coxph_model)
    # Extract p-values from zph_test
    p_values <- zph_test$table[3, "p"]
    return(p_values)
  }, error = function(e) {
    message("Error occurred: ", conditionMessage(e))
    return("NA")  # Return NULL or handle error case as needed
  })
}

pvalues_zph_test_circ_sex <- lapply(circ_sex_coxph, FUN = function(x) pvalue_zph_test_circ_sex(x))


names(pvalues_zph_test_circ_sex) <- circ_sex_coxph

## adjust pv

pvadj_zph_test_circ_sex <- p.adjust(as.numeric(unlist(pvalues_zph_test_circ_sex)), method="fdr")      
names(pvadj_zph_test_circ_sex) <- circ_sex_coxph            

coxmodel_result <- cbind(na.omit(pvadj_sex_coxph[pvadj_sex_coxph < 0.05]),pvadj_zph_test_circ_sex)
colnames(coxmodel_result) <- c("pvadj_coxph","pvadj_zphtest")

save_file_coxmodel_rest <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/07_survival_analysis/",project,"/coxph_zph_test_pvadj_sex_",project,".rds",sep="")

saveRDS(coxmodel_result, file=save_file_coxmodel_rest)                                      
                                    
}                                        
