library(dplyr)
library(tidyverse)
library(readr)
library(stringr)
library(readxl)

getwd()
setwd("C:/Users/lkgbe/OneDrive - Iowa State University/Documents/ISU/SPRING 2025/BCB4560/BCB456_Project/")
old_GeneID_S3 <- read_excel("./Table S3.xlsx")
head(old_GeneID_S3$geneID)
# Assuming your vector is old_GeneID$geneID
split_genes_S3 <- strsplit(old_GeneID_S3$geneID, "/")
all_geneIDs_S3 <- unique(unlist(split_genes_S3))
head(all_geneIDs_S3)
df_geneIDs_S3 <- data.frame(geneID = all_geneIDs_S3)
colnames(df_geneIDs_S3)[1] = "old"
write.table(df_geneIDs_S3, file = "all_geneIDs_S3.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


old_GeneID_S4 <- read_excel("./Table S4.xlsx")
head(old_GeneID_S4$geneID)
# Assuming your vector is old_GeneID$geneID
split_genes_S4 <- strsplit(old_GeneID_S4$geneID, "/")
all_geneIDs_S4 <- unique(unlist(split_genes_S4))
head(all_geneIDs_S4)
df_geneIDs_S4 <- data.frame(geneID = all_geneIDs_S4)
colnames(df_geneIDs_S4)[1] = "old"
write.table(df_geneIDs_S4, file = "all_geneIDs_S4.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
