args=(commandArgs(TRUE))


project <- args[1]
drug_target <- args[2]


library(gtools)
library(reshape2)


function_calc_ks_pvadj_logFC <- function(drug_target, project, sex){

	patients_geneexp_clind_drug_samptype <- readRDS(paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/02_selection_samples/02_1_ext_information_patients_tcga/",project,"/patients_geneexp_clind_drug_samptype_",project,".rds",sep=""))

pat_drug_db_circ <- readRDS(paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/04_mm_ko_drug_used/04_2_pat_drug_db_circ/",project,"/pat_drug_db_circ_",project,".rds",sep=""))

pathvals_drug <- readRDS(paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/04_mm_ko_drug_used/04_3_pathvals_ko_drug_used/",project,"/pathvals_",project,"_",drug_target,".rds",sep=""))

pathvals <- readRDS(paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/01_mechanistics_models/01_1_mechanistics_models_tcga/",project,"/pathvals_",project,".rds",sep=""))

drug_circs_belong <- pat_drug_db_circ[which(pat_drug_db_circ$entrez_id == drug_target),"circs"]

samples_sex <- patients_geneexp_clind_drug_samptype[patients_geneexp_clind_drug_samptype$gender == sex & patients_geneexp_clind_drug_samptype$sample_type == "Primary Solid Tumor","patient"]


## we apply ks.test and get the p.value
	ks_p.value <- lapply(1:length(drug_circs_belong), FUN=function(k) ks.test(pathvals_drug[drug_circs_belong[k],samples_sex],
			pathvals[drug_circs_belong[k],samples_sex])$p.value)
## we adjust pvalue by FDR	
	ks_p.adjust <- p.adjust(ks_p.value, method="fdr", n=length(ks_p.value))

	fc <- unlist(lapply(1:length(drug_circs_belong), FUN=function(k) mean(foldchange2logratio(foldchange(pathvals_drug[drug_circs_belong[k],samples_sex],
		pathvals[drug_circs_belong[k],samples_sex]), base=2))))

	df <- data.frame("kspvadj"=ks_p.adjust,"logFC"=fc, "drug_target"=rep(drug_target, length(drug_circs_belong)), "circ"=drug_circs_belong, "sex"=rep(sex, length(drug_circs_belong)),"project"=rep(project, length(drug_circs_belong)))
	
	return(df)
}

kstest_logfc_male <- function_calc_ks_pvadj_logFC(drug_target,project,"MALE")
kstest_logfc_female <- function_calc_ks_pvadj_logFC(drug_target,project,"FEMALE")

file_save <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/05_kstest_pvadj_logfc_kovsctrl/",project,"/kstest_logfc_",drug_target,".rds",sep="")

saveRDS(rbind(kstest_logfc_male, kstest_logfc_female), file=file_save)
