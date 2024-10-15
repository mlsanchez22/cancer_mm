#!/usr/bin/env Rscript

# Purpose: simplify drug actions, indicate frequent actions and add ATC descriptions.
# Input 1: path to drugbank parsed and translated table
# Input 2: path to descriptions of ATC categories
# Output: extended drugbank table with simplified actions and ATC descriptions

suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))

# get arguments from command line

args <- commandArgs(trailingOnly = TRUE)

# check if there are at least two arguments, otherwise return an error
if (length(args) < 2) {
  stop("Paths to drugbank and ATC tables required.", call. = F)
}

# load resources

message("Loading drugbank and ATC resources...")

drugbank_file <- args[1]
atc_file <- args[2]

drugbank_db <- read.delim(file = drugbank_file, sep = "\t" )
ATC_raw <- read.delim(atc_file, sep = ",") 

## check required columns

drugbank_cols <- c("name", "gene_symbol", "actions", "atc_codes")
if (!all(drugbank_cols %in% colnames(drugbank_db))) {
  stop(paste("Missing drugbank columns:", drugbank_cols[!(drugbank_cols %in% colnames(drugbank_db))]))
  }
atc_cols <- c("Class.ID", "ATC.LEVEL", "Preferred.Label")
if (!all(atc_cols %in% colnames(ATC_raw))) {
  stop(paste("Missing ATC columns:", atc_cols[!(atc_cols %in% colnames(ATC_raw))]))
}

# Reduce the number of drug actions

message("Simplifying drug actions...")

## unify the not-known drug actions to the same term "unknown"

drugbank_db <- drugbank_db %>% 
  mutate(actions = ifelse(str_detect(actions, "^$|other"), "unknown", actions))

##  Simplify the drug action vector to 5 categories

pos_eff <- c("activator", "(?>^\\b|\\||al )agonist", "inducer", "potentiator", 
             "stimulator", 
             "substrate", "ligand", 
             "allosteric")
neg_eff <- c("inhibit[^s]", "(?>se |ant)agonist", "blocker", "antisense", 
             "negat", "suppressor")
is_bind <- c("bind", "carrier", "chaperone")
is_ab <- "antibody"

actions_df <- drugbank_db %>% 
  group_by(actions) %>% summarize(n_drugs = n()) %>% ungroup() %>% 
  mutate(n_pos = unlist(lapply(actions, function(x) sum(str_detect(x, pos_eff)))), 
         n_neg = unlist(lapply(actions, function(x) sum(str_detect(x, neg_eff)))), 
         is_binder = unlist(lapply(actions, function(x) any(str_detect(x, is_bind)))), 
         is_antibody = unlist(lapply(actions, function(x) any(str_detect(x, is_ab))))) %>% 
  mutate(sum_ef = n_pos - n_neg) %>% 
  mutate(Drug_effect = case_when(sum_ef > 0 ~ "Activator", 
                                 sum_ef < 0 ~ "Inhibitor", 
                                 is_antibody == TRUE ~ "Antibody", 
                                 is_binder == TRUE ~ "Binder", 
                                 .default = "Other"), .before = 1)

## Add simplified effect to the database

drugbank_db <- drugbank_db %>% 
  left_join(dplyr::select(actions_df, actions, Drug_effect), by = "actions")


# aggregation of drug actions by gene or drug
## Since a drug can have several effects depending on the KDT it targets, and 
## KDTs can be targeted by many drugs, we will select the most common effect 
## that drugs have on each KDT to simplify plotting drug effects.

message("Aggregating drug actions...")

# aggregate drug actions by gene

drug_bygenes <- drugbank_db %>% 
  group_by(gene_symbol) %>% 
  summarize(Gene_effect_simpl = names(sort(table(Drug_effect), decreasing = T)[1]))

drugbank_db <- drugbank_db %>% left_join(drug_bygenes, by = "gene_symbol")


## Similarly, the various effects of drugs may be summarized
## if a drug has both activator and inhibitor effects, it will be denoted as Modulator

sing_eff <- drugbank_db %>% dplyr::select(name, Drug_effect) %>% distinct() %>% 
  group_by(name) %>% filter(n() == 1) %>%   # select drugs with a single drug effect
  mutate(Drug_effect_simpl = Drug_effect) %>% dplyr::select(name, Drug_effect_simpl) %>% unique()
multi_eff <- drugbank_db %>% dplyr::select(name, Drug_effect) %>% distinct() %>% 
  group_by(name) %>% filter(n() != 1, Drug_effect != "Other") %>%  # drugs with multiple effects and remove "Other"
  mutate(Drug_effect_simpl = case_when(n() == 1 ~ Drug_effect, 
                                       "Activator" %in% Drug_effect & "Inhibitor" %in% Drug_effect ~ "Modulator", 
                                       "Activator" %in% Drug_effect ~ "Activator", 
                                       "Inhibitor" %in% Drug_effect ~ "Inhibitor", 
                                       .default = "Modulator")) %>% 
  dplyr::select(name, Drug_effect_simpl) %>% unique()

drug_bydrugs <- rbind(sing_eff, multi_eff)

drugbank_db <- drugbank_db %>% left_join(drug_bydrugs, by = "name")



# HANDLE TABLE OF CATEGEORIES FROM ATC  
## Downloaded from https://bioportal.bioontology.org/ontologies/ATC

message("Incorporating ATC classification labels...")

ATC_raw$Class.ID <- str_split(ATC_raw$Class.ID, "ATC/") %>% sapply(., function(x) x[2])
ATC_raw <- ATC_raw[which(!is.na(ATC_raw$ATC.LEVEL)), ] %>% .[order(.$Class.ID),]

# fill blanks in drugbank ATC codes with drug synonyms

drugbank_db$synonyms =  sapply(str_split(drugbank_db$name, " ") , function(x) x[1])
noATC_all <- drugbank_db$atc_codes == ""
drugbank_db$atc_codes[noATC_all] <- drugbank_db$atc_codes[match(drugbank_db$synonyms[noATC_all],drugbank_db$synonyms[!noATC_all])]
drugbank_db$atc_codes[which(drugbank_db$atc_codes == "")] <- NA

drugbank_db <- dplyr::select(drugbank_db, -synonyms)

# explode ATC codes into one line each

drugbank_db <- drugbank_db %>% tidyr::separate_rows(atc_codes, sep = "\\|")

# create lists of substrings for the different ATC categories in drugbank

drugbank_db <- drugbank_db %>% mutate(
  atc_level_1_id = str_extract(atc_codes, "^."), 
  atc_level_1_label = ATC_raw$Preferred.Label[match(atc_level_1_id, ATC_raw$Class.ID)],
  atc_level_2_id = str_extract(atc_codes, "^.{3}"), 
  atc_level_2_label = ATC_raw$Preferred.Label[match(atc_level_2_id, ATC_raw$Class.ID)],
  atc_level_3_id = str_extract(atc_codes, "^.{4}"), 
  atc_level_3_label = ATC_raw$Preferred.Label[match(atc_level_3_id, ATC_raw$Class.ID)],
  atc_level_4_id = str_extract(atc_codes, "^.{5}"), 
  atc_level_4_label = ATC_raw$Preferred.Label[match(atc_level_4_id, ATC_raw$Class.ID)],
)


# write annotated table

message("Writing annotated drugbank table...")
drugbank_db <- drugbank_db %>% mutate_all(\(x) str_replace_all(x, "\n", " "))
write_tsv(drugbank_db, "drugbank_db-simpl_atc.tsv", quote = "needed")

message("Annotation of drugbank database done.")

