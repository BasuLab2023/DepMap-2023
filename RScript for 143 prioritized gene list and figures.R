##PIPELINE FOR 143 PROIRITIZED GENES##

#load and Install libraries
library(faraway)
library(MASS)
library(lars)
library(pls)
library(factoextra)
library(data.table)
library(dplyr)
library(readr)
library(tools)
library(purrr)
library(tidyverse)
library(magrittr)
library(sn)
library(pbapply)
library(mosaic)
library(matrixTests)
library(matrixStats)
library(plyr)
library(janitor)
library(MASS) # for chisq
library(descr) # for crosstable
library(ggpmisc)
library(ggrepel)
library(EnhancedVolcano)
library(lsr)
library(rstatix)
library(ggtext)
library(glue)
library(ggpubr)
library(grid)
library(formattable)
library(readxl)
library(gdata)
library(ggbreak)
library(ggplot2)


customGreen0 = "#DeF7E9"
customGreen = "#71CA97"
customRed = "#ff7f7f"

#download files
setwd("C:/Users/malay/Documents/Austin Data/Files")##Use the directory where the files are downloaded
gene_dependency_q2 <- read.csv("CRISPR_gene_dependency_23q2.csv", header = TRUE)
gene_effect_q2 <- read.csv("CRISPR_gene_effect_23q2.csv", header = TRUE)
gene_list_q2 <- read.csv("gene_list_23q2.csv", header = TRUE, check.names=FALSE)
gene_dependency <- read.csv("CRISPR_gene_dependency.csv", header = TRUE)
gene_effect <- read.csv("CRISPR_gene_effect.csv", header = TRUE)
HPV_neg_lines <- read.csv("HPV_neg_cell_lines.csv", header = TRUE)
DGI <- read.csv("categories_updated.csv", header = TRUE)
DGI_interaction <- read.csv("interactions.csv", header = TRUE)
CGC <- read.csv("CGC.csv", header = TRUE)
CGC_hnk <- read.csv("CGC_hnk.csv", header = TRUE)
core_fitness <- read.csv("core_fitness_genes.csv", header = TRUE)
CCLE.gene.cn <- fread("CCLE_gene_cn.csv")
CCLE.mutations <- fread("CCLE_mutations.csv")
HPVneg.cohort <- fread("HPVnegSCC-cohort.csv")
TCGA.mutations <- fread("Mutated_Genes_TCGA.txt")
TCGA.CNA <- fread("CNA_Genes_TCGA.txt")
oncokb_all <- read.csv("oncokb_genes.csv")
oncokb_not_onco <- oncokb_all %>% filter(Is.Oncogene == "No")
oncokb <- read.csv("oncokb_oncogenes.csv")
hallmark <- read.csv("hallmark.csv")
drugs <- read.csv("open_targets.csv")
drugs_window <- read.csv("open_targets_window.csv")
gene_names <- read.csv("all_gene_names.csv", header = TRUE)
HNSCC_gene_import <- read.csv("HNSCC_gene_import.csv", header = TRUE)
cancer_gene_import <- read.csv("cancer_gene_import.csv", header = TRUE)
hnscc_cell_lines <- read.csv("hnscc_cell_lines.csv", header = TRUE)
mutation_subgroups_cosmic <- read.csv("mutation_subgroups_cosmic.csv", header = TRUE)
cna_subgroups_cosmic <- read.csv("cna_subgroups_cosmic.csv", header = TRUE)
cna_cell_lines <- read.csv("cna_cell_lines.csv", header = TRUE)
cell_line_info <- read.csv("cell_line_info.csv", header = TRUE)
text_mining_gene_list <- read.csv("text_mining_gene_list.csv", header = TRUE)
text_mining_gene_list$genes <- str_trim(text_mining_gene_list$genes) 
CRISPR_common_essentials <- read.csv("CRISPR_common_essentials.csv")
david_groups <- read.csv("david_groups.csv")
cell_lines_recheck <- read.csv("cell_lines_recheck.csv")
hnscc_cell_lines_recheck <- read.csv("hnscc_cell_lines_recheck.csv")
original_list <- read.csv("original_list.csv")
hpv_pos_cell_lines <- read.csv("hpv_pos_cell_lines.csv")
open_targets_old <- read_excel("open_targets_old.xlsx")

#setup of cell lines and dependency data
HPV_neg_dependency <- setDT(gene_dependency)[DepMap_ID %chin% cell_lines_recheck$DepMap_ID]
gene_dependency_HNSCC <- HPV_neg_dependency %>% filter(DepMap_ID %in% hnscc_cell_lines_recheck$DepMap_ID)
gene_dependency_ESCC <- HPV_neg_dependency %>% filter(!(DepMap_ID %in% hnscc_cell_lines_recheck$DepMap_ID))
HPV_neg_dependency_cell_lines <- HPV_neg_dependency %>% select(DepMap_ID)
write.csv(HPV_neg_dependency_cell_lines, "all_cell_lines.csv")
HPV_neg_effect <- setDT(gene_effect)[DepMap_ID %chin% cell_lines_recheck$DepMap_ID]
HPV_neg_lines <- HPV_neg_dependency %>% select(DepMap_ID)
cell_lines_info <- cell_lines_recheck %>% filter(DepMap_ID %in% HPV_neg_lines$DepMap_ID)
write.csv(cell_lines_info, "cell_lines_info.csv")

HPV_neg_dependency_old <- setDT(gene_dependency)[DepMap_ID %chin% original_list$depmap_id]
HPV_neg_dependency_filter <- HPV_neg_dependency %>% filter(!(DepMap_ID %in% HPV_neg_dependency_old$DepMap_ID))
HPV_neg_dependency_old <- HPV_neg_dependency_old %>% filter(!(DepMap_ID %in% HPV_neg_dependency$DepMap_ID))

#setup of cell lines and dependency data (23q2)
names(gene_dependency_q2)[1] <- "DepMap_ID"
names(gene_effect_q2)[1] <- "DepMap_ID"
HPV_neg_dependency_q2 <- setDT(gene_dependency_q2)[DepMap_ID %chin% cell_lines_recheck$DepMap_ID]
gene_dependency_HNSCC_q2 <- HPV_neg_dependency_q2 %>% filter(DepMap_ID %in% hnscc_cell_lines_recheck$DepMap_ID)
gene_dependency_ESCC_q2 <- HPV_neg_dependency_q2 %>% filter(!(DepMap_ID %in% hnscc_cell_lines_recheck$DepMap_ID))
HPV_neg_dependency_cell_lines_q2 <- HPV_neg_dependency_q2 %>% select(DepMap_ID)
write.csv(HPV_neg_dependency_cell_lines_q2, "all_cell_lines_q2.csv")
HPV_neg_effect_q2 <- setDT(gene_effect_q2)[DepMap_ID %chin% cell_lines_recheck$DepMap_ID]
HPV_neg_lines_q2 <- HPV_neg_dependency_q2 %>% select(DepMap_ID)
cell_lines_info_q2 <- cell_lines_recheck %>% filter(DepMap_ID %in% HPV_neg_lines_q2$DepMap_ID)
write.csv(cell_lines_info_q2, "cell_lines_info_q2.csv")

HPV_neg_dependency_old_q2 <- setDT(gene_dependency_q2)[DepMap_ID %chin% original_list$depmap_id]
HPV_neg_dependency_filter_q2 <- HPV_neg_dependency_q2 %>% filter(!(DepMap_ID %in% HPV_neg_dependency_old_q2$DepMap_ID))
HPV_neg_dependency_old_q2 <- HPV_neg_dependency_old %>% filter(!(DepMap_ID %in% HPV_neg_dependency_q2$DepMap_ID))

#DATA PROCESSING------------------------------------------------------

#CNA data: Filter for cell lines in our cohort
CCLE.gene.cn.new <- CCLE.gene.cn %>% filter(V1 %in% HPV_neg_dependency$DepMap_ID)
write.csv(CCLE.gene.cn.new, "CCLE.cn.csv")

#Dependency data: Transpose so genes are row names
HPV_neg_dependency_t <- t(HPV_neg_dependency)
HPV_neg_effect_t <- t(HPV_neg_effect)
CCLE.gene.cn.new_t <- t(CCLE.gene.cn.new)

#Dependency data: Transpose so genes are row names (23q2)
HPV_neg_dependency_q2_t <- t(HPV_neg_dependency_q2)
HPV_neg_effect_q2_t <- t(HPV_neg_effect_q2)

#Dependency data: convert to numeric dataframe
HPV_neg_dependency_t.df <- data.frame(HPV_neg_dependency_t)
HPV_neg_effect_t.df <- data.frame(HPV_neg_effect_t)
CCLE.gene.cn.new_t.df <- data.frame(CCLE.gene.cn.new_t)

HPV_neg_dependency_t.n <- HPV_neg_dependency_t.df %>% mutate_all(as.numeric)
HPV_neg_effect_t.n <- HPV_neg_effect_t.df %>% mutate_all(as.numeric)
CCLE.gene.cn.new_t.n <- CCLE.gene.cn.new_t.df %>% mutate_all(as.numeric)

#Dependency data: convert to numeric dataframe (23q2)
HPV_neg_dependency_q2_t.df <- data.frame(HPV_neg_dependency_q2_t)
HPV_neg_effect_q2_t.df <- data.frame(HPV_neg_effect_q2_t)

HPV_neg_dependency_q2_t.n <- HPV_neg_dependency_q2_t.df %>% mutate_all(as.numeric)
HPV_neg_effect_q2_t.n <- HPV_neg_effect_q2_t.df %>% mutate_all(as.numeric)

#Dependency data: Filter if one essential cell line
HPV_neg_dependency_essential.n <- HPV_neg_dependency_t.n %>% filter_all(any_vars(. > 0.5))

HPV_neg_dependency_essential.n[HPV_neg_dependency_essential.n >= 0.5] <- 1
HPV_neg_dependency_essential.n[HPV_neg_dependency_essential.n < 0.5] <- 0
HPV_neg_dependency_essential.n[is.na(HPV_neg_dependency_essential.n)] <- 0
HPV_neg_dependency_essential.n <- HPV_neg_dependency_essential.n %>% mutate_all(as.numeric)
HPV_neg_dependency_essential.n$dependent_cell_lines <- rowSums(HPV_neg_dependency_essential.n)

HPV_neg_dependency_essential.n <- HPV_neg_dependency_essential.n %>% filter(dependent_cell_lines > 0)
HPV_neg_dependency_essential.n <- HPV_neg_dependency_essential.n[,-c(88)]

#Dependency data: Filter if one essential cell line (23q2)
HPV_neg_dependency_q2_essential.n <- HPV_neg_dependency_q2_t.n %>% filter_all(any_vars(. > 0.5))

HPV_neg_dependency_q2_essential.n[HPV_neg_dependency_q2_essential.n >= 0.5] <- 1
HPV_neg_dependency_q2_essential.n[HPV_neg_dependency_q2_essential.n < 0.5] <- 0
HPV_neg_dependency_q2_essential.n[is.na(HPV_neg_dependency_q2_essential.n)] <- 0
HPV_neg_dependency_q2_essential.n <- HPV_neg_dependency_q2_essential.n %>% mutate_all(as.numeric)
HPV_neg_dependency_q2_essential.n$dependent_cell_lines <- rowSums(HPV_neg_dependency_q2_essential.n)

HPV_neg_dependency_q2_essential.n <- HPV_neg_dependency_q2_essential.n %>% filter(dependent_cell_lines > 0)
HPV_neg_dependency_q2_essential.n <- HPV_neg_dependency_q2_essential.n[,-c(88)]

#Dependency data: cell lines as column names
names <- HPV_neg_effect$DepMap_ID
colnames(HPV_neg_dependency_t.n) <- names
colnames(HPV_neg_effect_t.n) <- names

colnames(HPV_neg_dependency_essential.n) <- names
write.csv(HPV_neg_dependency_essential.n, "HPV_neg_dependency_essential.csv", row.names = TRUE)

names_cn <- CCLE.gene.cn.new$V1
colnames(CCLE.gene.cn.new_t.n) <- names_cn

#Dependency data: cell lines as column names (23q2)
names <- HPV_neg_effect_q2$DepMap_ID
colnames(HPV_neg_dependency_q2_t.n) <- names
colnames(HPV_neg_effect_q2_t.n) <- names

colnames(HPV_neg_dependency_q2_essential.n) <- names
write.csv(HPV_neg_dependency_q2_essential.n, "HPV_neg_dependency_essential.csv", row.names = TRUE)

#clean up gene names in cn table
CCLE.cn.1 <- rownames(CCLE.gene.cn.new_t.n)
CCLE.cn.2 <- gsub("\\s*\\([^\\)]+\\)","",as.character(CCLE.cn.1))
CCLE.cn.2 <- data.frame(CCLE.cn.2)
colnames(CCLE.cn.2) <- "genes"
CCLE_cn_nice <- cbind(CCLE.cn.2, CCLE.gene.cn.new_t.n)

#filter for cn levels
CCLE_cn_levels <- CCLE_cn_nice
genes <- CCLE_cn_levels[,1]
CCLE_cn_levels[,1] <- NULL

CCLE_cn_levels <- CCLE_cn_levels %>% mutate_all(~ case_when(. >= 1.214 ~ "amp",
                                                            . <= 0.731 ~ "del",
                                                            TRUE ~ 'none'))

CCLE_cn_levels_nice <- cbind(genes, CCLE_cn_levels)

#Dependency data (value): clean up gene names
HPV_neg_dependency_essential_2 <- rownames(HPV_neg_dependency_t.n)
HPV_neg_dependency_essential_3 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_essential_2)
HPV_neg_dependency_essential_4 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_essential_3)
HPV_neg_dependency_rows <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_essential_4)
write.table(HPV_neg_dependency_rows, file = "HPV_neg_dependency_rows.csv", sep = "\t",
            row.names = FALSE, quote = FALSE) 

HPV_neg_dependency_rows.df <- data.frame(HPV_neg_dependency_rows)
colnames(HPV_neg_dependency_rows.df) <- "genes"
HPV_neg_dependency_nice <- cbind(HPV_neg_dependency_rows.df, HPV_neg_dependency_t.n)

#Dependency data (value): clean up gene names (23q2)
HPV_neg_dependency_q2_essential_2 <- rownames(HPV_neg_dependency_q2_t.n)
HPV_neg_dependency_q2_essential_3 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_q2_essential_2)
HPV_neg_dependency_q2_essential_4 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_q2_essential_3)
HPV_neg_dependency_q2_rows <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_q2_essential_4)
write.table(HPV_neg_dependency_q2_rows, file = "HPV_neg_dependency_q2_rows.csv", sep = "\t",
            row.names = FALSE, quote = FALSE) 

HPV_neg_dependency_q2_rows.df <- data.frame(HPV_neg_dependency_q2_rows)
colnames(HPV_neg_dependency_q2_rows.df) <- "genes"
HPV_neg_dependency_nice_q2 <- cbind(HPV_neg_dependency_q2_rows.df, HPV_neg_dependency_q2_t.n)

#Dependency data (binary): clean up gene names, export essential genes
HPV_neg_dependency_essential_2 <- rownames(HPV_neg_dependency_essential.n)
HPV_neg_dependency_essential_3 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_essential_2)
HPV_neg_dependency_essential_4 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_essential_3)
HPV_neg_dependency_essential_rows <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_essential_4)
write.table(HPV_neg_dependency_essential_rows, file = "HPV_neg_dependency_essential_rows", sep = "\t",
            row.names = FALSE, quote = FALSE) 

HPV_neg_dependency_essential_rows.df <- data.frame(HPV_neg_dependency_essential_rows)
colnames(HPV_neg_dependency_essential_rows.df) <- "genes"
HPV_neg_dependency_essential_nice <- cbind(HPV_neg_dependency_essential_rows.df, HPV_neg_dependency_essential.n)

HPV_neg_dependency_2 <- rownames(HPV_neg_dependency_t.n)
HPV_neg_dependency_3 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_2)
HPV_neg_dependency_4 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_3)
HPV_neg_dependency_rows <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_4)
write.table(HPV_neg_dependency_rows, file = "HPV_neg_dependency_rows", sep = "\t",
            row.names = FALSE, quote = FALSE) 

#Dependency data (binary): clean up gene names, export essential genes (23q2)
HPV_neg_dependency_q2_essential_2 <- rownames(HPV_neg_dependency_q2_essential.n)
HPV_neg_dependency_q2_essential_3 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_q2_essential_2)
HPV_neg_dependency_q2_essential_4 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_q2_essential_3)
HPV_neg_dependency_q2_essential_rows <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_q2_essential_4)
write.table(HPV_neg_dependency_q2_essential_rows, file = "HPV_neg_dependency_q2_essential_rows", sep = "\t",
            row.names = FALSE, quote = FALSE) 

HPV_neg_dependency_q2_essential_rows.df <- data.frame(HPV_neg_dependency_q2_essential_rows)
colnames(HPV_neg_dependency_q2_essential_rows.df) <- "genes"
HPV_neg_dependency_essential_nice_q2 <- cbind(HPV_neg_dependency_q2_essential_rows.df, HPV_neg_dependency_q2_essential.n)

HPV_neg_dependency_q2_2 <- rownames(HPV_neg_dependency_q2_t.n)
HPV_neg_dependency_q2_3 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_q2_2)
HPV_neg_dependency_q2_4 <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_q2_3)
HPV_neg_dependency_q2_rows <- sub("^(.*)[.].*", "\\1", HPV_neg_dependency_q2_4)
write.table(HPV_neg_dependency_q2_rows, file = "HPV_neg_dependency_q2_rows", sep = "\t",
            row.names = FALSE, quote = FALSE) 

#Dependency data: Filter for druggable genome
DGI_druggable <- dplyr::filter(DGI, grepl('DRUGGABLE GENOME',category))
HPV_neg_druggable <- setDT(HPV_neg_dependency_essential_nice)[genes %chin% DGI_druggable$entrez_gene_symbol]

write.table(HPV_neg_druggable$genes, file = "HPV_neg_druggable", sep = "\t",
            row.names = FALSE, quote = FALSE) 

#Dependency data: Filter for druggable genome (23q2)
HPV_neg_druggable_q2 <- setDT(HPV_neg_dependency_essential_nice_q2)[genes %chin% DGI_druggable$entrez_gene_symbol]

write.table(HPV_neg_druggable_q2$genes, file = "HPV_neg_druggable_q2", sep = "\t",
            row.names = FALSE, quote = FALSE) 

#RELATIVE DEPENDENCY vs GENE EFFECT ----------------------------------------
#Reset to 23q2 dataset (activate to run)
HPV_neg_druggable <- HPV_neg_druggable_q2
HPV_neg_effect_t.n <- HPV_neg_effect_q2_t.n
HPV_neg_dependency_nice <- HPV_neg_dependency_nice_q2

#Dependency data:: mark gene dependency > 0.5 in cohort only, calculated relative dependency
HPV_neg_druggable_combined_mark <- HPV_neg_druggable
HPV_neg_druggable_genes <- HPV_neg_druggable[,1]
HPV_neg_druggable_combined_mark[HPV_neg_druggable_combined_mark >= 0.5] <- 1
HPV_neg_druggable_combined_mark[HPV_neg_druggable_combined_mark < 0.5] <- 0
HPV_neg_druggable_combined_mark[,1] <- NULL
HPV_neg_druggable_combined_mark[is.na(HPV_neg_druggable_combined_mark)] <- 0
HPV_neg_druggable_combined_mark.n <- HPV_neg_druggable_combined_mark %>% mutate_all(as.numeric)

HPV_neg_druggable_combined_mark_nice <- cbind(HPV_neg_druggable_genes, HPV_neg_druggable_combined_mark)
HPV_neg_druggable_combined_mark_nice$dependent_cell_lines <- rowSums(HPV_neg_druggable_combined_mark.n)
HPV_neg_druggable_combined_mark_nice$relative_dependency <- (HPV_neg_druggable_combined_mark_nice$dependent_cell_lines)/(ncol(HPV_neg_druggable_combined_mark_nice)-2)

write.csv(HPV_neg_druggable_combined_mark_nice, file = "HPV_neg_druggable_combined_mark_nice.csv",
          row.names = TRUE) 

#Check core fitness against OncoKB and druggable
core_onco <- setDT(core_fitness)[genes %chin% oncokb$Gene]
core_not_onco <- setDT(core_fitness)[!(genes %chin% oncokb$Gene)]
core_onco_druggable <- setDT(HPV_neg_druggable_combined_mark_nice)[genes %chin% core_onco$genes]
core_druggable <- setDT(HPV_neg_druggable_combined_mark_nice)[genes %chin% core_fitness$genes]

#Dependency data: filter core fitness genes
HPV_neg_no_core <- setDT(HPV_neg_druggable_combined_mark_nice)[!(genes %chin% core_fitness$genes)]

#filter common essentials
CRISPR_common_essentials <- as.data.frame(lapply(CRISPR_common_essentials,function(x) if(is.character(x)|is.factor(x)) gsub(" .*","",x) else x))
HPV_neg_no_common <- setDT(HPV_neg_druggable_combined_mark_nice)[!(genes %chin% CRISPR_common_essentials$gene)]
HPV_neg_no_core_common <- setDT(HPV_neg_no_core)[!(genes %chin% CRISPR_common_essentials$gene)]

#Dependency data: filter relative dependency > 10% in cohort
hist(HPV_neg_no_core_common$dependent_cell_lines, breaks=76,main="Number of cell lines with gene dependency")
favstats(HPV_neg_no_core_common$dependent_cell_lines)

HPV_neg_druggable_combined_mark_nice_RD <- setDT(HPV_neg_no_core_common)[dependent_cell_lines > 7]

write.csv(HPV_neg_druggable_combined_mark_nice_RD, file = "HPV_neg_druggable_combined_mark_nice_RD.csv",
          row.names = TRUE) 

#Dependency data: filter RD <= 9%
RD_less_than_9 <- setDT(HPV_neg_no_core_common)[dependent_cell_lines < 8 & dependent_cell_lines > 1]
write.csv(RD_less_than_9, file = "RD_less_than_9.csv",
          row.names = TRUE) 
hist(RD_less_than_9$dependent_cell_lines)

#clean up gene names in HPV neg gene effect table
HPV_neg_effect_essential_2 <- rownames(HPV_neg_effect_t.n)
HPV_neg_effect_essential_3 <- sub("^(.*)[.].*", "\\1", HPV_neg_effect_essential_2)
HPV_neg_effect_essential_4 <- sub("^(.*)[.].*", "\\1", HPV_neg_effect_essential_3)
HPV_neg_effect_rows <- sub("^(.*)[.].*", "\\1", HPV_neg_effect_essential_4)
write.table(HPV_neg_effect_rows, file = "HPV_neg_effect_rows", sep = "\t",
            row.names = FALSE, quote = FALSE) 

HPV_neg_effect_rows.df <- data.frame(HPV_neg_effect_rows)
colnames(HPV_neg_effect_rows.df) <- "genes"
HPV_neg_effect_nice <- cbind(HPV_neg_effect_rows.df, HPV_neg_effect_t.n)

#median gene effect in essential cell lines, in only gene list
HPV_neg_druggable_effect <- setDT(HPV_neg_effect_nice)[genes %chin% HPV_neg_druggable_combined_mark_nice_RD$genes]
Effect_dependent_only <- HPV_neg_druggable_effect
Effect_dependent_only[,c(2:88)][Effect_dependent_only[,c(2:88)] > -0.5] <- NA
Effect_dependent_only.n <- data.matrix(Effect_dependent_only)
Effect_dependent_only$median_gene_effect <- rowMedians(Effect_dependent_only.n[,c(2:88)], na.rm=TRUE)

#median gene effect in essential cell lines, in druggable
HPV_neg_druggable_effect2 <- setDT(HPV_neg_effect_nice)[genes %chin% HPV_neg_druggable$genes]
Effect_dependent_only2 <- HPV_neg_druggable_effect2

Effect_dependent_only2_correct <- HPV_neg_druggable_effect2[,c(2:88)]*HPV_neg_druggable_combined_mark_nice[,c(2:88)]
Effect_dependent_only2_correct_nice <- cbind(HPV_neg_druggable_effect2[,c(1)],Effect_dependent_only2_correct)
Effect_dependent_only2_correct_nice[,c(2:88)][Effect_dependent_only2_correct_nice[,c(2:88)] == 0] <- NA
Effect_dependent_only2_correct_nice.n <- data.matrix(Effect_dependent_only2_correct_nice)
Effect_dependent_median_gene_effect <- as.data.frame(rowMedians(Effect_dependent_only2_correct_nice.n[,c(2:88)], na.rm=TRUE))
colnames(Effect_dependent_median_gene_effect) <- "median_gene_effect"
Effect_dependent_only2_done <- cbind(Effect_dependent_only2_correct_nice,Effect_dependent_median_gene_effect)

Effect_dependent_only2[,c(2:88)][Effect_dependent_only2[,c(2:88)] > -0.5] <- NA
Effect_dependent_only2.n <- data.matrix(Effect_dependent_only2)
Effect_dependent_only2$median_gene_effect <- rowMedians(Effect_dependent_only2.n[,c(2:88)], na.rm=TRUE)

#Invert effect scores
HPV_neg_druggable_effect[,c(2:88)] <- -1*HPV_neg_druggable_effect[,c(2:88)]
HPV_neg_druggable_effect2[,c(2:88)] <- -1*HPV_neg_druggable_effect2[,c(2:88)]
Effect_dependent_only[,c(2:89)] <- -1*Effect_dependent_only[,c(2:89)]
Effect_dependent_only2[,c(2:89)] <- -1*Effect_dependent_only2[,c(2:89)]
Effect_dependent_only2_done[,c(2:89)] <- -1*Effect_dependent_only2_done[,c(2:89)]

#print gene probabilities (binary and values) in gene list
HPV_neg_druggable_dependency <- setDT(HPV_neg_dependency_nice)[genes %chin% HPV_neg_druggable_combined_mark_nice_RD$genes]

HPV_neg_druggable_dependency$median = rowMedians(as.matrix(HPV_neg_druggable_dependency[,c(2:88)], na.rm=TRUE))

write.csv(HPV_neg_druggable_dependency,"HPV_neg_druggable_dependency.csv")

HPV_neg_druggable_dependency_binary <- setDT(HPV_neg_dependency_essential_nice)[genes %chin% HPV_neg_druggable_combined_mark_nice_RD$genes]
write.csv(HPV_neg_druggable_dependency_binary,"HPV_neg_druggable_dependency_binary.csv")

#median gene probability in essential cell lines
Probability_dependent_only <- HPV_neg_druggable
Probability_dependent_only[,c(2:88)][Probability_dependent_only[,c(2:88)] < 0.5] <- NA
Probability_dependent_only.n <- data.matrix(Probability_dependent_only)
Probability_dependent_only$median <- rowMedians(Probability_dependent_only.n[,c(2:88)], na.rm=TRUE)

#plot RD vs. gene effect in all druggable genes
RDplot <- merge(HPV_neg_druggable_combined_mark_nice, Effect_dependent_only2_done, by="genes")
RDplot$relative_dependency = RDplot$relative_dependency * 100

RDplot$common <- RDplot$genes %in% unlist(CRISPR_common_essentials$gene) %>% as.character() %>% str_to_title()
RDplot$core <- RDplot$genes %in% unlist(core_fitness$genes) %>% as.character() %>% str_to_title()
RDplot$list <- RDplot$genes %in% unlist(HPV_neg_druggable_combined_mark_nice_RD$genes) %>% as.character() %>% str_to_title()

RDplot <- RDplot %>% mutate(common_category = case_when(common == 'True' & core == 'True' ~ 'Common essential and core fitness  (n = 143)',
                                                        common == 'True'  ~ 'Common essential gene only (n = 98)',
                                                        core == 'True'  ~ 'Core fitness gene only (n = 59)',
                                                        TRUE ~ 'Not common essential or core fitness (n = 720)'))
table(RDplot$common_category)

RDplot <- RDplot %>% mutate(category = case_when(common == 'True' & core == 'True' ~ 'Common essential and core fitness  (n = 143)',
                                                 common == 'True'  ~ 'Common essential gene only (n = 98)',
                                                 core == 'True'  ~ 'Core fitness gene only (n = 59)',
                                                 list == 'False' & core == 'False' ~ 'Essential in < 10% cell lines (n = 581)',
                                                 list == 'True' ~ 'Essential in > 10% cell lines (n = 139)'))
table(RDplot$category)

RDplot <- RDplot %>% mutate(category_3 = case_when(common == 'True' | core == 'True' ~ 'Common essential or core fitness (n = 303)',
                                                   list == 'False' & core == 'False' ~ 'Essential in < 8 cell lines (n = 704)',
                                                   list == 'True' ~ 'Prioritized druggable targets (n = 143)'))
table(RDplot$category_3)

RDplot <- RDplot %>% mutate(core_binary = case_when(core == 'True' ~ 'Core fitness genes (n = 202)',
                                                    TRUE ~ 'Not core fitness genes (n = 818)'
))
table(RDplot$core_binary)

RDplot <- RDplot %>% mutate(rare_binary = case_when(list == 'False' & core == 'False' ~ 'Essential in < 10% cell lines (n = 581)',
                                                    list == 'True' ~ 'Essential in > 10% cell lines (n = 139)'
))
table(RDplot$rare_binary)
RDplot$rare_binary <- factor(RDplot$rare_binary, levels=c('Essential in < 10% cell lines (n = 581)','Essential in > 10% cell lines (n = 139)'))

#RDplot <- RDplot %>% mutate(Cancer = case_when(HNSCC == 'True' ~ 'Clinical inhibitor in HNSCC', cancer == 'True' ~ 'Clinical inhibitor in cancer'))

RDplot_druggable <- RDplot 

RDplot_window <- subset(RDplot, median_gene_effect > 0.57)
RDplot_cancer <- setDT(RDplot)[genes %chin% cancer_gene_import$genes]
RDplot_inhibitor <- setDT(RDplot)[genes %chin% drugs$'Gene Symbol']
RDplot_HNSCC <- setDT(RDplot)[genes %chin% HNSCC_gene_import$genes]

RDplot$cancer <- RDplot$genes %in% unlist(RDplot_cancer$genes) %>% as.character() %>% str_to_title()
RDplot$HNSCC <- RDplot$genes %in% unlist(RDplot_HNSCC$genes) %>% as.character() %>% str_to_title()
table(RDplot_druggable$category_3)

#calculate number of common essentials
RDplot_strong <- subset(RDplot, median_gene_effect > 1)
RDplot_common <- subset(RDplot, relative_dependency > 90)
RDplot_common_strong <- subset(RDplot, relative_dependency > 90 & median_gene_effect > 1.005)
RDplot_common_strong_core_overlap <- subset(RDplot_common_strong, core_binary == 'Core fitness genes (n = 202)')
RDplot_common_core_overlap <- subset(RDplot_common, core_binary == 'Core fitness genes (n = 202)')

RDplot_common_strong_core_combo <- subset(RDplot, (relative_dependency > 90 & median_gene_effect > 1.005) | core_binary == 'Core fitness genes (n = 202)')

#quick check gene effect scores of common vs. core
RDplot_core_true <- RDplot %>% filter(core == "True")
RDplot_common_true <- RDplot %>% filter(common == "True")
RDplot_common_or_core_true <- RDplot %>% filter(common == "True" | core == "True")
RDplot_common_and_core_true <- RDplot %>% filter(common == "True" & core == "True")

RDplot_core_true %>% dplyr::summarise(median = median(median_gene_effect, na.rm=TRUE), iqr = quantile(median_gene_effect, c(0.25, 0.5, 0.75), na.rm=TRUE))
RDplot_common_true %>% dplyr::summarise(median = median(median_gene_effect, na.rm=TRUE), iqr = quantile(median_gene_effect, c(0.25, 0.5, 0.75), na.rm=TRUE))
RDplot_common_or_core_true %>% dplyr::summarise(median = median(median_gene_effect, na.rm=TRUE), iqr = quantile(median_gene_effect, c(0.25, 0.5, 0.75), na.rm=TRUE))

#####GENERATE FIGURE 1B

RDplot$median_gene_effect <- -1*(RDplot$median_gene_effect)

text <- paste("Essential in <10% cell lines")
fig2a_pre <- ggplot(RDplot, aes(x=relative_dependency, y=median_gene_effect)) + 
  scale_y_reverse() +
  annotate("rect",xmin=-Inf,xmax=9,ymin=-Inf,ymax=Inf,alpha=.1,color="gold",fill="gold") +
  geom_point(shape=19, size=1.7, alpha=0.7, aes(colour=category_3)) + 
  labs(title="", x="% Cell lines with essentiality", y="Median gene effect score") +
  #geom_smooth(method = "lm", formula = my.formula, color="grey50", se=FALSE, size=0.8) +
  theme_classic() +
  theme(legend.title = element_blank(), text = element_text(size = 15)) +
  scale_fill_manual(drop = FALSE) +
  scale_color_manual(values=c("firebrick1","dodgerblue","chartreuse3")) 
print(fig2a_pre)

write.csv(RDplot, file = "RDplot.csv", row.names = TRUE)

#calculate spearmans for 2A
cor.test(x=RDplot$relative_dependency, y=RDplot$median_gene_effect, method = 'spearman')

#Figure 2B: highlight strong/rare
RDplot <- RDplot %>% filter(core == 'False' & common == 'False')
table(RDplot$rare_binary)
strong_rare_RDplot <- subset(RDplot, median_gene_effect > 0.57 & rare_binary == 'Essential in 3-10% cell lines (n = 155)')

#how many oncogenes in lower?
RDplot_onco <- setDT(RDplot)[genes %chin% oncokb$Gene]
RDplot_onco_rare <- RDplot_onco %>% filter(rare_binary == "Essential in 3-9% cell lines (n = 155)")

fig2b_pre <- ggplot(RDplot, aes(x=relative_dependency, y=median_gene_effect)) + 
  geom_point(shape=1, size=1.7, alpha=0.7, aes(colour=rare_binary)) + 
  labs(x="Percent cell lines", y="Median gene effect score") +
  geom_text_repel(data = subset(RDplot, genes == ''), na.rm = TRUE, hjust = -1, force = 10, nudge_x = 5, nudge_y = .05, segment.color = 'grey', box.padding = 0.2, point.padding = 0.2) +
  theme_bw() +
  theme(legend.title = element_blank(), text = element_text(size = 15)) +
  scale_color_manual(values=c("firebrick1","dodgerblue"))

print(fig2b_pre)
fig2b <- ggplotGrob(fig2b_pre)
fig2b$layout$clip[fig2b$layout$name=="panel"] <- "off"
grid.draw(fig2b)
write.csv(RDplot, file = "RDplot.csv", row.names = TRUE)

#categories box plot
RDplot_druggable %>% group_by(core) %>% dplyr::summarise(median = median(median_gene_effect, na.rm=TRUE), iqr = quantile(median_gene_effect, c(0.25, 0.5, 0.75), na.rm=TRUE))
#RDplot_not_core_druggable  %>% group_by(rare) %>% dplyr::summarise(median = median(median_gene_effect, na.rm=TRUE), iqr = quantile(median_gene_effect, c(0.25, 0.5, 0.75), na.rm=TRUE))
RDplot_common %>% dplyr::summarise(median = median(median_gene_effect, na.rm=TRUE), iqr = quantile(median_gene_effect, c(0.25, 0.5, 0.75), na.rm=TRUE))

RDplot_druggable %>% group_by(category_3) %>% dplyr::summarise(median = median(median_gene_effect, na.rm=TRUE), iqr = quantile(median_gene_effect, c(0.25, 0.5, 0.75), na.rm=TRUE))
RDplot_druggable %>% group_by(category_3) %>% shapiro_test(median_gene_effect)
KW.category_3 <- RDplot_druggable %>% kruskal_test(median_gene_effect ~ category_3)
KW.category_3
KWpair.category_3 <- RDplot_druggable %>% dunn_test(median_gene_effect ~ category_3, p.adjust.method = "bonferroni")
KWpair.category_3
formattable(KWpair.category_3, align ="c")

KWpair.category_3 <- KWpair.category_3 %>% add_xy_position(x = "category_3")
RDplot_druggable$median_gene_effect <- -1*(RDplot_druggable$median_gene_effect)

#FIGURE 1C FINAL
RDplot_druggable$category_3 <- factor(RDplot_druggable$category_3, levels = c('Common essential or core fitness (n = 303)', 'Essential in < 10% cell lines (n = 704)', 'Prioritized druggable targets (n = 143)'))
fig2c <- ggboxplot(RDplot_druggable, x = "category_3", y = "median_gene_effect", outlier.shape = NA) +
  stat_pvalue_manual(KWpair.category_3, hide.ns = FALSE, y.position = -3, tip.length = - 0.03, step.increase = 0.1) +
  geom_jitter(width = 0.20, shape=19, size=1.7, alpha=0.5, aes(colour=category_3)) + 
  labs(title = "Median Gene Effect Scores of Druggable Dependencies", x = "Categories of genes", y = "Median gene effect score", subtitle = get_test_label(KW.category_3, detailed = TRUE),
       caption = get_pwc_label(KWpair.category_3)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2), labels = NULL) +
  theme(axis.title.x = element_blank(), legend.position="none", text = element_text(size = 15)) +
  scale_color_manual(values=c("firebrick1","dodgerblue","chartreuse3"))
ggpar(fig2c, ylim = c(-0.3,-3.5))

#plot RD vs. gene effect in gene list
RDplot2 <- merge(HPV_neg_druggable_combined_mark_nice_RD, Effect_dependent_only2_done, by="genes")
RDplot2$relative_dependency = RDplot2$relative_dependency * 100
RDplot2$core <- RDplot2$genes %in% unlist(core_fitness$genes) %>% as.character() %>% str_to_title()
RDplot2$list <- RDplot2$genes %in% unlist(HPV_neg_druggable_combined_mark_nice_RD$genes) %>% as.character() %>% str_to_title()
RDplot2 <- RDplot2 %>% mutate(Category = case_when(list == 'True' ~ 'Genes in > 10% cell lines and not core fitness',
                                                   core == 'True' ~ 'Core fitness genes',
                                                   TRUE ~ 'Genes in < 10% cell lines'))
RDplot2 <- RDplot2 %>% mutate(Regression = case_when(relative_dependency > 99 ~ 'Genes in 100% cell lines',
                                                     TRUE ~ 'Genes in 10-99% cell lines'))

#filter to new gene list
RDplot2 <- RDplot2 %>% filter(relative_dependency > 3 & core == 'False') #& !(relative_dependency > 90 & median_gene_effect > 1.005))

my.formula <- y ~ poly(x, 2)
form <- lm(median_gene_effect ~ relative_dependency,data=RDplot2)
dat <- predict(form, interval="confidence", level= 0.99)
#RDplot2$inside <- ifelse(RDplot2$median_gene_effect < dat[,"upr"] & RDplot2$median_gene_effect > dat[,"lwr"], "", as.character(RDplot2$genes))
RDplot2$inside <- ifelse(RDplot2$median_gene_effect > 0.57, "", as.character(RDplot2$genes))

RDplot2_cancer <- setDT(RDplot2)[genes %chin% cancer_gene_import$genes]
RDplot2_inhibitor <- setDT(RDplot2)[genes %chin% drugs$'Gene Symbol']

RDplot2_HNSCC <- setDT(RDplot2)[genes %chin% HNSCC_gene_import$genes]
RDplot2_window <- subset(RDplot2, median_gene_effect > 0.57)
RDplot2_window <- RDplot2_window %>% filter(!(genes == 'EGFR' | genes == 'TP63' | genes == 'TYMS' | genes == 'ERBB3' | genes == 'PIK3CA' | genes == 'ERBB2' | genes == 'CDK6' | genes == 'MCL1' | genes == 'AURKA'))
HNSCC <- subset(RDplot2, genes == 'TUBB4B' | genes == 'EGFR' | genes == 'CDK6' | genes == 'TYMS' | genes == 'ERBB3' | genes == 'PIK3CA' | genes == 'ERBB2' | genes == 'TP63' | genes == 'MCL1' | genes == 'AURKA' | genes == 'AURKA' | genes == 'AURKA')
#TUBB4B, TUBA1C, TUBA1B, TUBB, NDUFA4, NDUFS3
options(ggrepel.max.overlaps = Inf)

RDplot2$cancer <- RDplot2$genes %in% unlist(RDplot2_cancer$genes) %>% as.character() %>% str_to_title()
RDplot2$HNSCC <- RDplot2$genes %in% unlist(RDplot2_HNSCC$genes) %>% as.character() %>% str_to_title()
RDplot2 <- RDplot2 %>% mutate(Cancer = case_when(HNSCC == 'True' ~ 'Clinical inhibitor in HNSCC',
                                                 cancer == 'True' ~ 'Clinical inhibitor in cancer'))


RDplot2 <- RDplot2 %>% mutate(dplot = case_when(median_gene_effect > 0.57 ~ 'Median gene effect score >= EGFR (n = 24)',
                                                TRUE ~ 'Median gene effect score < EGFR (n = 115)'))
table(RDplot2$dplot)

histogram(RDplot2$median_gene_effect, breaks = 20)

#calculate spearmans for 2D
cor.test(x=RDplot2$relative_dependency, y=RDplot2$median_gene_effect, method = 'spearman')

table(RDplot2$HNSCC)
RDplot2 %>% group_by(HNSCC) %>% dplyr::summarise(median = median(median_gene_effect, na.rm=TRUE), iqr = quantile(median_gene_effect, c(0.25, 0.5, 0.75)))
RDplot2 %>% group_by(HNSCC) %>% shapiro_test(median_gene_effect)
KW.HNSCC <- RDplot2 %>% kruskal_test(median_gene_effect ~ HNSCC)
KW.HNSCC
KWpair.HNSCC <- RDplot2 %>% dunn_test(median_gene_effect ~ HNSCC, p.adjust.method = "bonferroni")
KWpair.HNSCC

RDplot2 %>%
  group_by(HNSCC) %>%
  ggplot(., aes(x=HNSCC, y=`median_gene_effect`)) + 
  geom_boxplot(width=0.3) + 
  geom_point(alpha = 2/10) + 
  ylim(.5,2.0) + 
  labs(title = "Median Gene Effect Scores", x = "Categories of genes", y = "Median gene effect score")

#GENE LIST------------------------------------------------------
#open targets
drugs <- read.csv("open_targets.csv")
colnames(drugs) <- c("Gene Symbol", "FDA Approved Indications in Cancer","FDA Approved Indications", "Phase II Indications", "Inhibitor ≥ Phase II",	"Bucket 1: FDA Approved Inhibitor?", "FDA Approved Inhibitors",	"Bucket 2: Phase II/III Inhibitors?",	"Phase II/III Inhibitors", "Bucket 3: Phase I or Preclinical Inhibitor?",	 "Clinical Inhibitor in Cancer",	"HNSCC Inhibitor Phase II or Better?", "Clinical Inhibitor in HNSCC")
drugs <- drugs[!apply(drugs == "", 1, all),]

#gene list
gene_list <- data.frame(cbind(HPV_neg_druggable_combined_mark_nice$genes, HPV_neg_druggable_combined_mark_nice$dependent_cell_lines, HPV_neg_druggable_combined_mark_nice$relative_dependency, Effect_dependent_only2_done$median_gene_effect))
colnames(gene_list) <- c("genes","dependent_cell_lines","relative_dependency","median_gene_effect")
gene_list <- arrange(gene_list,desc(relative_dependency),desc(median_gene_effect))


HPV_neg_drug_interaction <- setDT(HPV_neg_druggable_combined_mark_nice_RD)[genes %chin% DGI_interaction$gene_name]
DGI_interaction_gene_list <- setDT(DGI_interaction)[gene_name %chin% gene_list$genes]
DGI_interaction_gene_list <- DGI_interaction_gene_list[,c("gene_name","interaction_types","drug_name","drug_claim_primary_name","interaction_group_score")] %>% distinct(gene_name,drug_name, .keep_all = TRUE) %>% arrange(gene_name,desc(interaction_group_score))
colnames(DGI_interaction_gene_list)[colnames(DGI_interaction_gene_list) == 'gene_name'] <- 'genes'
DGI_interaction_gene_list <- DGI_interaction_gene_list %>% filter(interaction_types == "inhibitor" | interaction_types == "antagonist" | interaction_types == "modulator")
DGI_interaction_gene_list_top_score <- DGI_interaction_gene_list[!duplicated(DGI_interaction_gene_list$genes),]
gene_list$drug_interaction <- gene_list$genes %in% unlist(HPV_neg_drug_interaction$genes) %>% as.character() %>% str_to_title()

DGI_names <- DGI[,c("entrez_gene_symbol", "gene_long_name")] %>% distinct()
colnames(DGI_names) <- c("genes", "gene_long_name")

DGI_MF <- DGI %>% filter(!category == 'DRUGGABLE GENOME' & !category == 'CLINICALLY ACTIONABLE' & !category == 'DRUG RESISTANCE')
DGI_MF$category <- DGI_MF$category %>% fct_relevel("SERINE THREONINE KINASE","TYROSINE KINASE","PROTEASE","PROTEASE INHIBITOR","DNA REPAIR","TRANSCRIPTION FACTOR","ION CHANNEL","G PROTEIN COUPLED RECEPTOR","NUCLEAR HORMONE RECEPTOR","HISTONE MODIFICATION","ABC TRANSPORTER") 
DGI_MF <- DGI_MF[,c("entrez_gene_symbol","category")] %>% arrange(entrez_gene_symbol,category)
colnames(DGI_MF)[colnames(DGI_MF) == 'entrez_gene_symbol'] <- 'genes'
DGI_MF_top <- DGI_MF[!duplicated(DGI_MF$genes),]

gene_list <- dplyr::inner_join(gene_list, DGI_names, by="genes")
gene_list <- dplyr::left_join(gene_list, DGI_interaction_gene_list_top_score, by="genes")
gene_list <- dplyr::left_join(gene_list, DGI_MF_top, by="genes")
gene_list$known_cancer_oncokb <- gene_list$genes %in% unlist(oncokb$Gene) %>% as.character() %>% str_to_title()

#finish putting together gene list
gene_list <- arrange(gene_list,desc(relative_dependency),desc(median_gene_effect))
gene_list$gene_long_name <- str_to_sentence(gene_list$gene_long_name)
gene_list$category <- str_to_sentence(gene_list$category)
gene_list$interaction_types <- str_to_sentence(gene_list$interaction_types)
gene_list$drug_claim_primary_name <- str_to_sentence(gene_list$drug_claim_primary_name)

gene_list$relative_dependency = digits(gene_list$relative_dependency, digits=0)
gene_list$relative_dependency = gene_list$relative_dependency * 100

gene_list$median_gene_effect = digits(gene_list$median_gene_effect, digits=2)
gene_list[is.na(gene_list)] <- ""

gene_names_simple <- gene_names[c("genes","names")]
gene_list <- left_join(gene_list, gene_names_simple, by="genes")
gene_list <- left_join(gene_list, text_mining_gene_list, by="genes")

gene_list_all_columns <- gene_list
gene_list <- gene_list[,c("relative_dependency","median_gene_effect","genes", "names","known_cancer_oncokb", "interaction_types","Cancer","Oncogene","HNSCC")]
gene_list <- gene_list %>% mutate(known_cancer_oncokb = str_replace(known_cancer_oncokb, "True", "Yes"))
gene_list <- gene_list %>% mutate(known_cancer_oncokb = str_replace(known_cancer_oncokb, "False", "No"))

gene_list <- gene_list %>% mutate(interaction_types = str_replace(interaction_types, "Inhibitor", "Yes"))
gene_list <- gene_list %>% mutate(interaction_types = str_replace(interaction_types, "Antagonist", "Yes"))
gene_list <- gene_list %>% mutate(interaction_types = str_replace(interaction_types, "Modulator", "Yes"))
gene_list$interaction_types[gene_list$interaction_types==""] <- "No" 

gene_list$clinical_inhibitor <- ifelse(gene_list$genes %in% drugs$`Gene Symbol`, 'Yes', 'No')

#subset the gene list and finalize tables
gene_list$median_gene_effect = digits(gene_list$median_gene_effect, digits=2)
gene_list_druggable <- gene_list
gene_list <- setDT(gene_list)[genes %chin% RDplot2$genes]
gene_list <- gene_list %>% distinct()
strong_list <- setDT(gene_list)[`median_gene_effect` > 1]

#strong_rare special table
gene_list_all_columns <- gene_list_all_columns[,c("dependent_cell_lines","median_gene_effect","genes", "names","known_cancer_oncokb", "interaction_types","Cancer","Oncogene","HNSCC")]
gene_list_all_columns <- gene_list_all_columns %>% mutate(known_cancer_oncokb = str_replace(known_cancer_oncokb, "True", "Yes"))
gene_list_all_columns <- gene_list_all_columns %>% mutate(known_cancer_oncokb = str_replace(known_cancer_oncokb, "False", "No"))

gene_list_all_columns <- gene_list_all_columns %>% mutate(interaction_types = str_replace(interaction_types, "Inhibitor", "Yes"))
gene_list_all_columns <- gene_list_all_columns %>% mutate(interaction_types = str_replace(interaction_types, "Antagonist", "Yes"))
gene_list_all_columns <- gene_list_all_columns %>% mutate(interaction_types = str_replace(interaction_types, "Modulator", "Yes"))
gene_list_all_columns$interaction_types[gene_list_all_columns$interaction_types==""] <- "No" 

gene_list_all_columns$clinical_inhibitor <- ifelse(gene_list_all_columns$genes %in% drugs$`Gene Symbol`, 'Yes', 'No')

gene_list_all_columns$median_gene_effect = digits(gene_list_all_columns$median_gene_effect, digits=2)
gene_list_all_columns_druggable <- gene_list_all_columns
gene_list_all_columns <- setDT(gene_list_all_columns)[genes %chin% HPV_neg_druggable_combined_mark_nice_RD$genes]

strong_rare_list <- setDT(as.data.frame(gene_list_all_columns_druggable))[genes %chin% strong_rare_RDplot$genes]

core_onco_list <- setDT(as.data.frame(gene_list_druggable))[genes %chin% core_onco_druggable$genes]
window_list <- gene_list %>% filter(median_gene_effect > 0.57)

#column names and reorder
colnames(gene_list) <- c("% Cell Lines","Median Gene Effect","Gene Symbol", "Gene Name","Known Oncogene","Known Inhibitor","Cancer","Oncogene","HNSCC","Clinical Inhibitor")
colnames(strong_rare_list) <- c("# Cell Lines","Median Gene Effect","Gene Symbol", "Gene Name","Known Oncogene","Known Inhibitor","Cancer","Oncogene","HNSCC", "Clinical Inhibitor")
colnames(core_onco_list) <- c("% Cell Lines","Median Gene Effect","Gene Symbol", "Gene Name","Known Oncogene","Known Inhibitor","Cancer","Oncogene","HNSCC", "Clinical Inhibitor")
colnames(window_list) <- c("% Cell Lines","Median Gene Effect","Gene Symbol", "Gene Name","Known Oncogene","Known Inhibitor","Cancer","Oncogene","HNSCC", "Clinical Inhibitor")

gene_list <- gene_list[,c("% Cell Lines","Median Gene Effect","Gene Symbol", "Gene Name","Known Oncogene","Known Inhibitor","Clinical Inhibitor")]
strong_rare_list <- strong_rare_list[,c("# Cell Lines","Median Gene Effect","Gene Symbol", "Gene Name","Known Oncogene","Known Inhibitor","Clinical Inhibitor")]
core_onco_list <- core_onco_list[,c("% Cell Lines","Median Gene Effect","Gene Symbol", "Gene Name","Known Oncogene","Known Inhibitor","Clinical Inhibitor")]
window_list <- window_list[,c("% Cell Lines","Median Gene Effect","Gene Symbol", "Gene Name","Known Oncogene","Known Inhibitor","Clinical Inhibitor")]
final_window_list <- window_list %>% filter(`Known Oncogene` == 'No')
known_window_list <- window_list %>% filter(`Known Oncogene` == 'Yes')

#fix hist names for gene list
gene_list$`Gene Symbol` <- str_replace(gene_list$`Gene Symbol`, "HIST1H1C", "H1-2")
gene_list$`Gene Symbol` <- str_replace(gene_list$`Gene Symbol`, "HIST2H3D", "H3C13")
gene_list$`Gene Symbol` <- str_replace(gene_list$`Gene Symbol`, "H3F3B", "H3-3B")
gene_list$`Gene Symbol` <- str_replace(gene_list$`Gene Symbol`, "H3F3A", "H3-3A")
gene_list$`Gene Symbol` <- str_replace(gene_list$`Gene Symbol`, "HIST1H3A", "H3C1")
gene_list$`Gene Symbol` <- str_replace(gene_list$`Gene Symbol`, "HIST1H3D", "H3C4")
gene_list$`Gene Symbol` <- str_replace(gene_list$`Gene Symbol`, "HIST1H3I", "H3C11")
gene_list$`Gene Symbol` <- str_replace(gene_list$`Gene Symbol`, "HIST1H4E", "H4C5")

#make formattables
gene_list$'Median Gene Effect' <- -1*(gene_list$'Median Gene Effect')
formattable(gene_list, align ="c")
formattable(core_onco_list, align ="c")
formattable(window_list, align ="c")
formattable(final_window_list, align ="c")
formattable(known_window_list, align ="c")

#update the known inhibitor list
gene_list_drugs <- gene_list %>% filter(`Known Inhibitor` == "Yes")
gene_list_drugs_old <- filter(gene_list_drugs, `Gene Symbol` %in% open_targets_old$`Gene`)
gene_list_drugs_new <- filter(gene_list_drugs, !(`Gene Symbol` %in% open_targets_old$`Gene`))
write.csv(gene_list_drugs_new, "gene_list_drugs_new.csv")

#gene_list_drugs_new <- setDT(gene_list_drugs)[!(`Gene Symbol` %chin% drugs$`Gene Symbol`)]
drugs_delete <- filter(drugs, !(`Gene Symbol` %in% gene_list$`Gene Symbol` ))

yike <- c("EGLN1",
          "P2RY6",
          "PCSK9",
          "UGCG",
          "BIRC2",
          "ITGB1",
          "MAP3K11")
drugs_old <- filter(gene_list_drugs, `Gene Symbol` %in% yike)
write.csv(drugs_old, "drugs_old.csv")

open_targets_done <- setDT(open_targets_old)[`Gene` %chin% gene_list_drugs$`Gene Symbol`]

#update the window list
window_genes_curation <- read_excel("window genes curation.xlsx")
window_genes_curation$`Gene Symbol`
window_list_old <- setDT(window_list)[`Gene Symbol` %chin% window_genes_curation$`Gene Symbol`]
window_list_new <- setDT(window_list)[!(`Gene Symbol` %chin% window_genes_curation$`Gene Symbol`)]
window_list_delete <- setDT(window_genes_curation)[!(`Gene Symbol` %chin% window_list$`Gene Symbol`)]
window_list_drugs <- window_list %>% filter(`Known Inhibitor` == "Yes")
window_list_drugs_old <- setDT(window_list)[`Gene Symbol` %chin% drugs$`Gene Symbol`]
window_list_drugs_new_1 <- setDT(window_list)[!(`Gene Symbol` %chin% drugs$`Gene Symbol`)]
window_list_drugs_new <- window_list_drugs_new_1 %>% filter(`Known Inhibitor` == "Yes")

write.csv(window_list_drugs, file = "window_list_drugs.csv", row.names = TRUE)
write.csv(window_list_old, file = "window_list_old.csv", row.names = TRUE)

write.csv(window_list_new, file = "window_list_new.csv", row.names = TRUE)
write.csv(known_window_list, file = "known_window_list.csv", row.names = TRUE)
write.csv(final_window_list, file = "final_window_list.csv", row.names = TRUE)
write.csv(window_list, file = "window_list.csv", row.names = TRUE)
write.csv(strong_rare_list, file = "strong_rare_list.csv", row.names = TRUE)
write.csv(gene_list, file = "gene_list.csv", row.names = TRUE)
write.csv(gene_list_druggable, file = "gene_list_druggable.csv", row.names = TRUE)
write.csv(core_onco_list, file = "core_onco_list.csv", row.names = TRUE)

####OUTPUT IN THE MANUSCRIPT###

###FINAL TABLE OUTPUT as in SUPPLEMENTARY TABLE 1###

formattable(gene_list, align ="c")

#################END##########################################################################
