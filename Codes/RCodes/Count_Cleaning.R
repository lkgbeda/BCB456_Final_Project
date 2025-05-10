# Load necessary library
library(tidyverse)

# Set input and output directories
input_dir <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Merged/HTSeq_Count"
output_dir <- file.path(input_dir, "/Cleaned_Count")



# List all HTSeq count files
#sample_files <- list.files(input_dir, pattern = "SRR\\d+_counts\\.txt", full.names = TRUE)
sample_files <- list.files(input_dir, pattern = "merged_(CD|ND)[0-9]+_counts\\.txt", full.names = TRUE)

# Process each file
for (file in sample_files) {
  # Read the file
  count_data <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  
  # Extract the 2nd and 3rd columns (V2 and V3)
  cleaned_data <- count_data[, c(1, 3)]
  
  # Get the sample name
  #sample_name <- stringr::str_extract(basename(file), "SRR\\d+")
  sample_name <- stringr::str_extract(basename(file), "(CD|ND)[0-9]+")
  
  # Define output file path
  output_file <- file.path(output_dir, paste0(sample_name, "_counts.txt"))
  
  # Write the cleaned data to the new file
  write.table(cleaned_data, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}
