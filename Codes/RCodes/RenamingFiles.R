# Load required library
library(readr)
library(dplyr)
library(stringr)

# Paths
csv_path <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/metadata.csv"
input_dir <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Coordinate_Sorted"
output_dir <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Unmerged"

# Ensure output directory exists
#if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Read CSV
metadata <- read_csv(csv_path)

# List all files in the input directory
files <- list.files(input_dir, full.names = TRUE)

# Iterate through each file
for (file_path in files) {
  file_name <- basename(file_path)
  
  # Extract SRR ID prefix from filename using regex
  srr_id <- str_extract(file_name, "SRR[0-9]+")
  
  # Match it to the metadata
  match_row <- metadata %>% filter(`Accession - SRR` == srr_id)
  
  if (nrow(match_row) == 1) {
    # Construct new file name
    group <- match_row$Group
    suffix <- str_remove(file_name, srr_id)  # get the remaining part like '.bam' or '.bam.bai'
    new_name <- paste0(srr_id, "_", group, suffix)
    new_path <- file.path(output_dir, new_name)
    
    # Copy and rename the file
    file.copy(file_path, new_path)
    cat("Renamed:", file_name, "->", new_name, "\n")
  } else {
    cat("No match found for:", file_name, "\n")
  }
}
