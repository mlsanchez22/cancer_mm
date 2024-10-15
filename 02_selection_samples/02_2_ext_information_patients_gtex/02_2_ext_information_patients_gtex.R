### selecction of patients with clinical data sex and age
### creating a subset of gene expression

# Patient information and extrated place sample
annotation_sampleattributes <- read.table("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", sep="\t", header=T)

# dataset with gene expresion
gtex_counts <- read.table("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct", skip=2, head=T)
colnames(gtex_counts) <- gsub(".","-", colnames(gtex_counts), fixed=T)

annotation_sampleattributes_comp <- annotation_sampleattributes[annotation_sampleattributes$SAMPID %in% colnames(gtex_counts),]
annotation_sampleattributes_comp <- annotation_sampleattributes_comp[annotation_sampleattributes_comp$SMTSD != "",]

# patient information with sex and age
annotation_subjectphenotypes <- read.table("/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt", sep="\t", header=T)
rownames(annotation_subjectphenotypes) <- annotation_subjectphenotypes$SUBJID

annotation_sampleatributes_sel <- annotation_sampleattributes[,c("SAMPID","SMTS","SMTSD")]
rownames(annotation_sampleatributes_sel) <- annotation_sampleatributes_sel$SAMPID


# extract id of patients
subjib_patients_comp <- paste("GTEX-",unique(unlist(lapply(strsplit(annotation_sampleattributes_comp$SAMPID,"-", fixed=T), FUN=function(x) x[2]))), sep="")

# creation of dataframe with patients and annotations about sex and age

annotation_subjectphenotypes_comp <- annotation_subjectphenotypes[subjib_patients_comp,]
list_patient <- lapply(annotation_subjectphenotypes_comp$SUBJID, FUN=function(x) annotation_sampleatributes_sel$SAMPID[grep(x, annotation_sampleatributes_sel$SAMPID)])
list_sex <- do.call(list,lapply(1:length(annotation_subjectphenotypes_comp$SEX), function(x) rep(annotation_subjectphenotypes_comp$SEX[x],length(unlist(list_patient[x])))))    
list_age <- do.call(list,lapply(1:length(annotation_subjectphenotypes_comp$AGE), function(x) rep(annotation_subjectphenotypes_comp$AGE[x],length(unlist(list_patient[x])))))
df_patient_sex_age <- data.frame(patients = unlist(list_patient),
                                   sex_b = unlist(list_sex),
                                   age = unlist(list_age))
rownames(df_patient_sex_age) <- df_patient_sex_age$patients

df_patient_sex_age_comp <- df_patient_sex_age[annotation_sampleattributes_comp$SAMPID,]
df_patient_sex_age_comp$sex <- rep("-", length(df_patient_sex_age_comp$patients))
df_patient_sex_age_comp$sex[which(df_patient_sex_age_comp$sex_b == 2)] <- "female"
df_patient_sex_age_comp$sex[which(df_patient_sex_age_comp$sex_b == 1)] <- "male"
df_patient_sex_age_comp$tissue <- annotation_sampleatributes_sel[df_patient_sex_age_comp$patients,"SMTS"]
annotation_comp <- df_patient_sex_age_comp[,c("patients","age","sex","tissue")]
annotation_comp <- annotation_comp[annotation_comp$tissue != "",]

saveRDS(annotation_comp, file="/mnt/lustre/scratch/CBRA/research/projects/tcga_mm/results/02_selection_samples/02_2_ext_information_patients_gtex/donors_sex_age_tissue.rds")

