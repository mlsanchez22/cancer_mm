args=(commandArgs(TRUE))

drug_target <- args[1]

library(gtools)
library(reshape2)

function_calc_ks_pvadj_logFC <- function(drug_target, tissue, sex){

## we know with samples belong to each tissue and sex
	donors_sex_age_tissue <- readRDS("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/02_selection_samples/02_2_ext_information_patients_gtex/donors_sex_age_tissue.rds")

## we know with that belong to cicuits are in hipathia
	circ_physi_genes <- readRDS("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/04_mm_ko_drug_used/circ_physi_genes.rds")

## circuit activation of ko
	pathvals_drug <- readRDS(paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/04_mm_ko_drug_used/04_4_pathvals_ko_drug_used_gtex/pathvals_",drug_target,"_gtex.rds",sep=""))


## circuit activation of control
	pathvals <- readRDS("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/01_mechanistics_models/01_2_mechanistics_models_gtex/path_vals_gtex.rds")

## samples with sex and tissue
	samples <- donors_sex_age_tissue[intersect(which(donors_sex_age_tissue$tissue == tissue),which(donors_sex_age_tissue$sex == sex)),"patients"]

## circuits that belong gene
	circ_belong_gene <- circ_physi_genes[which(circ_physi_genes$genes == drug_target),"circs"]

## we apply ks.test and get the p.value
	ks_p.value <- lapply(1:length(circ_belong_gene), FUN=function(k) ks.test(pathvals_drug[circ_belong_gene[k],samples], pathvals[circ_belong_gene[k],samples])$p.value)
## we adjust pvalue by FDR	
	ks_p.adjust <- p.adjust(ks_p.value, method="fdr", n=length(ks_p.value))
	
	fc <- unlist(lapply(1:length(circ_belong_gene), FUN=function(k) mean(foldchange2logratio(foldchange(pathvals_drug[circ_belong_gene[k],samples], pathvals[circ_belong_gene[k],samples]), base=2))))
	df <- data.frame("kspvadj"=ks_p.adjust,"logFC"=fc, "ko"=rep(drug_target, length(circ_belong_gene)), "circ"=circ_belong_gene, "sex"=rep(sex, length(circ_belong_gene)),"tissue"=rep(tissue, length(tissue)))
	return(df)

}

donors_sex_age_tissue <- readRDS("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/02_selection_samples/02_2_ext_information_patients_gtex/donors_sex_age_tissue.rds")

circ_physi_genes <- readRDS("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/04_mm_ko_drug_used/circ_physi_genes.rds")

tissues_female <- unique(donors_sex_age_tissue[which(donors_sex_age_tissue$sex == "female"),"tissue"])
tissues_male <- unique(donors_sex_age_tissue[which(donors_sex_age_tissue$sex == "male"),"tissue"])

tissues <- setdiff(intersect(tissues_female, tissues_male), "Kidney")

results_female <- do.call(rbind,lapply(tissues, FUN=function(tissue) function_calc_ks_pvadj_logFC(drug_target, tissue, "female")))
results_male <- do.call(rbind,lapply(tissues, FUN=function(tissue) function_calc_ks_pvadj_logFC(drug_target, tissue, "male")))

results <- rbind(results_female, results_male)

file_save <- paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/05_kstest_pvadj_logfc_kovsctrl/gtex/kstest_pvadj_logfc_kovsctrl_",drug_target,".rds", sep="")
saveRDS(results, file= file_save)

