args=(commandArgs(TRUE))

tissue_ctrol <- args[1]

library(dplyr)

gtex_counts <- read.table("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2, head=T)
colnames(gtex_counts) <- gsub(".","-", colnames(gtex_counts), fixed=T)

donors_sex_age_tissue <- readRDS("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/02_selection_samples/02_2_ext_information_patients_gtex/donors_sex_age_tissue.rds")

gtex_counts_female <- gtex_counts[,as.data.frame(donors_sex_age_tissue %>% filter(tissue==tissue_ctrol & sex=="female"))$patients]
file_drugrep_gtex_female <- cbind(gtex_counts[,c("Name","Description")], gtex_counts_female)
write.table(file_drugrep_gtex_female, file = paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/09_drug_repository/09_1_genexpresion_gtex_tiss_ctrol/gtex_geneexp_female_",tissue_ctrol,".tsv",sep=""), row.names=FALSE, sep="\t",quote=F)

gtex_counts_male <- gtex_counts[,as.data.frame(donors_sex_age_tissue %>% filter(tissue==tissue_ctrol & sex=="male"))$patients]
file_drugrep_gtex_male <- cbind(gtex_counts[,c("Name","Description")], gtex_counts_male)
write.table(file_drugrep_gtex_male, file = paste("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/09_drug_repository/09_1_genexpresion_gtex_tiss_ctrol/gtex_geneexp_male_",tolower(tissue_ctrol),".tsv",sep=""), row.names=FALSE, sep="\t",quote=F)

