library(tidyverse)

# Paths
input_dir <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/HTSeq_Count"
output_dir <- file.path(input_dir, "GeneNaming")
mapping_file <- file.path(output_dir, "gene_id_to_name_mapping.tsv")

# Create output directory if needed
if (!dir.exists(output_dir)) dir.create(output_dir)

# List count files
sample_files <- list.files(input_dir, pattern = "SRR\\d+_counts\\.txt", full.names = TRUE)

# Initialize list for gene name mapping
gene_name_mapping <- list()

# Process each file
for (file in sample_files) {
  df <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  colnames(df) <- c("gene_id", "external_name", "count")
  
  # Keep only gene_id and count for DESeq2
  cleaned_df <- df[, c("gene_id", "count")]
  
  # Store gene mapping from one representative file
  if (length(gene_name_mapping) == 0) {
    gene_name_mapping <- df[, c("gene_id", "external_name")] %>% unique()
    
    # Modify external_name to keep only the second part (after comma)
    gene_name_mapping <- gene_name_mapping %>%
      mutate(external_name = sapply(strsplit(external_name, ","), function(x) trimws(x[2])))
  }
  
  # Output cleaned file
  sample_name <- str_extract(basename(file), "SRR\\d+")
  out_file <- file.path(output_dir, paste0(sample_name, "_cleaned.txt"))
  write.table(cleaned_df, out_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Save gene name mapping to a file
write.table(gene_name_mapping, mapping_file, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
