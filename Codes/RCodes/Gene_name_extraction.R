
library(dplyr)
library(tidyverse)
library(readr)
library(stringr)
library(readxl)

getwd()
setwd("C:/Users/lkgbe/OneDrive - Iowa State University/Documents/ISU/SPRING 2025/BCB4560/BCB456_Project/")

gff_data <- read_tsv("genomic.gff", comment = "#", col_names = FALSE)
head(gff_data)
tail(gff_data)

colnames(gff_data) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
head(gff_data)
tail(gff_data)
library(dplyr)
gff_data <- gff_data %>%
  mutate(gene_id = sub(".*ID=([^;]+);?.*", "\\1", attributes))
head(gff_data)
tail(gff_data)

# Extract specific tags from the attributes column
gff_data_extracted <- gff_data %>%
  mutate(
    locus_tag = str_extract(attributes, "locus_tag=[^;]+"),
    old_locus_tag = str_extract(attributes, "old_locus_tag=[^;]+"),
    Name = str_extract(attributes, "Name=[^;]+")
  ) %>%
  mutate(
    locus_tag = str_remove(locus_tag, "locus_tag="),
    old_locus_tag = str_remove(old_locus_tag, "old_locus_tag="),
    Name = str_remove(Name, "Name=")
  )
head(gff_data_extracted)
tail(gff_data_extracted)

geneifo <- gff_data_extracted[c("locus_tag","Name")]
colnames(geneifo)[1] <- "old"
head(geneifo)
print(geneifo %>%filter(is.na(locus_tag)))
new_tag <- read_excel("./Extra files/old to new names C58.xlsx")
head(new_tag)
new_tag <-new_tag%>%select(Aliases,old)
colnames(new_tag)[1] <- "GeneID"

# Check for duplicated values in each
duplicated_new_tag <- new_tag$old[duplicated(new_tag$old)]
duplicated_geneifo <- geneifo$old[duplicated(geneifo$old)] 

# Remove duplicates based on "old" in both data frames
new_tag_unique <- new_tag %>% distinct(old, .keep_all = TRUE)
geneifo_unique <- geneifo %>% distinct(old, .keep_all = TRUE)

new_geneinfo <- right_join(new_tag_unique, geneifo_unique[, c("old", "Name")], by = "old")
sum(is.na(new_geneinfo)[3])
write.table(new_geneinfo, file = "new_geneinfo.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
