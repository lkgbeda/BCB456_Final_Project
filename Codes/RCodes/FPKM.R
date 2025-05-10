
# Load required libraries for data processing and visualization
library(rtracklayer)  # For importing GTF files
library(dplyr)        # For data manipulation
library(tidyr)        # For reshaping data
library(ggplot2)      # For creating plots
library(pheatmap)     # For heatmap generation
library(readr)        # For reading delimited files

#-----------------------------#
# Set Paths and Directories
#-----------------------------#
# Define the base directory where your project files are located
base_dir <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project"

# Define paths to input and output files
gtf_file <- file.path(base_dir, "Ref_Genomes/Agrobacterium_fabrum/protein_coding_genes.gtf") # GTF file for annotation
count_dir <- file.path(base_dir, "HTSeq_Count/Cleaned_Count_Updated") # Folder containing HTSeq count files
fig_dir <- file.path(base_dir, "Figures") # Directory for saving plots
deg_dir <- file.path(base_dir, "DEA_DESeq") # Not used here but typically for differential expression results
metadata_file <- file.path(base_dir, "metadata.csv") # CSV file with sample metadata
gene_info_file <- file.path(base_dir, "new_geneinfo.tsv") # Gene name and ID mapping
fpkm_output <- file.path(base_dir, "HTSeq_Count/fpkm_values.txt") # Output file for FPKM values

#-----------------------------#
#Extract Gene Lengths from GTF
#-----------------------------#

# Read GTF
gtf <- import(gtf_file)

# For each gene, calculate gene length (you can use exon ranges summed if needed)
gene_lengths <- as.data.frame(gtf)


# In simple bacterial GTFs like Agrobacterium, one entry per gene might exist, so you can just:
gene_lengths_simple <- gene_lengths[, c("gene_id", "start", "end")]

# Calculate length
gene_lengths_simple$Length <- abs(gene_lengths_simple$end - gene_lengths_simple$start) + 1

# Keep only unique gene_id and length
gene_lengths_final <- gene_lengths_simple[, c("gene_id", "Length")]
colnames(gene_lengths_final) <- c("GeneID", "Length")

# Convert to kilobases
gene_lengths_final$Length_KB <- gene_lengths_final$Length / 1000

#----------------------------------#
# Load and Combine HTSeq Count Files
#----------------------------------#

# List all HTSeq count files
count_files <- list.files(count_dir, pattern = "SRR\\d+_counts\\.txt", full.names = TRUE)
# Load gene lengths (from the GTF file used in HTSeq, this can be extracted in a separate step)
# You can load gene lengths into a data frame 'gene_lengths' with columns for gene_id and length in base pairs.

 # Now read them in, renaming 'ReadCount' to the sample name (derived from filename)
count_data <- lapply(count_files, function(file) {
    counts <- read.table(file, header = FALSE, stringsAsFactors = FALSE)
    sample_name <- gsub("_count.txt", "", basename(file))  # Extract sample name from file name
    colnames(counts) <- c("GeneID", sample_name)  # Rename columns uniquely
    return(counts)
  })
  
# Merge all the tables by GeneID
combined_counts <- Reduce(function(x, y) merge(x, y, by = "GeneID", all = TRUE), count_data)
  
#----------------------------------#
# Calculate FPKM
#----------------------------------#

# Calculate total reads for each sample
total_reads <- colSums(combined_counts[, -1])  # Excluding GeneID column

str(gene_lengths)           # Check structure of the data frame
head(gene_lengths)          # Peek at the first few rows
sum(is.na(gene_lengths$Length))  # See if Length column exists and has NAs

# Assuming gene_lengths is a data frame with columns GeneID and Length (in base pairs)
gene_lengths$Length_KB <- gene_lengths$width / 1000 # Convert to kilobases

# Merge the gene lengths with the counts
combined_counts <- merge(combined_counts, gene_lengths_final, by = "GeneID", all.x = TRUE)

# Create a copy of combined_counts
fpkm <- combined_counts

# Find sample columns
sample_cols <- setdiff(colnames(fpkm), c("GeneID", "Length", "Length_KB"))

# Calculate FPKM
for (sample in sample_cols) {
  fpkm[[sample]] <- (fpkm[[sample]] / fpkm$Length_KB) / (total_reads[sample] / 1e6)
}

# Resulting fpkm matrix
head(fpkm)
tail(fpkm)
str(fpkm)
view(fpkm)
# Save the FPKM data frame to a text file (tab-separated, no row names or quotes)
write.table(fpkm, fpkm_output, sep = "\t", row.names = FALSE, quote = FALSE)

# Read the saved FPKM file back into R
fpkm <- read_table("/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/HTSeq_Count/fpkm_values.txt")

#----------------------------------#
# Load Metadata and Gene Info
#----------------------------------#
# Read the gene ID to gene name mapping file
gene_name_mapping <- read.delim(gene_info_file, sep = "\t", stringsAsFactors = FALSE)

# List all HTSeq count files that match the SRR pattern
sample_files <- list.files(count_dir, pattern = "SRR\\d+_counts\\.txt", full.names = FALSE)

# Preview the first and last few sample file names
head(sample_files)
tail(sample_files)

# Extract just the SRR sample IDs from the file names
sample_names <- stringr::str_extract(basename(sample_files), "SRR\\d+")
head(sample_names)
tail(sample_names)

# Read the metadata file that contains sample info (e.g., SRX/SRR, group, day)
sample_table <- read.table(metadata_file,
                           header = TRUE, sep = ",", fill = TRUE, quote = "\"",
                           comment.char = "", stringsAsFactors = FALSE)

# Add a column for the count file names to the sample table
sample_table$file_name <- sample_files

# Rename relevant columns for consistency
colnames(sample_table)[c(1,4)] <- c("SRX", "SRR")
head(sample_table)

# Reorder and filter sample table columns
sample_table <- sample_table[, c(2,8,7,1,6,5,2,3,4,6)] 

# Sort sample table by replicate
sample_table <- sample_table[order(sample_table$Replicate),]
head(sample_table)
tail(sample_table)

# Show the distribution of samples across groups and days
table(sample_table$Day, sample_table$Group)

# Convert Day and Group columns to factors
sample_table$Day <- factor(sample_table$Day)
sample_table$Group <- factor(sample_table$Group)

# View levels of Day and Group
levels(sample_table$Day)
levels(sample_table$Group)

# Ensure the gene ID column is named "GeneID" for joining
colnames(gene_name_mapping)[1] <- "GeneID"
head(gene_name_mapping)

#----------------------------------#
# Prepare FPKM for Plotting
#----------------------------------#

# Remove 'Length' and 'Length_KB' columns from FPKM for tidy formatting
fpkm_clean <- fpkm %>%
  dplyr::select(-Length, -Length_KB)

# Pivot FPKM data into long format (GeneID, file_name, FPKM)
fpkm_long <- pivot_longer(fpkm_clean, 
                          cols = -GeneID, 
                          names_to = "file_name", 
                          values_to = "FPKM")

# Preview the reshaped FPKM data
head(fpkm_long)
tail(fpkm_long)

# Add sample name and group info to the long FPKM table
fpkm_long <- fpkm_long %>%
  left_join(sample_table %>% dplyr::select(file_name, Sample_Name, Group), by = "file_name")

# Join gene names based on GeneID; fallback to GeneID if name is missing
fpkm_long <- fpkm_long %>%
  left_join(gene_name_mapping, by = "GeneID") %>%
  mutate(GeneName = ifelse(is.na(Name), GeneID, Name))

#----------------------------------#
# Plot FPKM Distributions
#----------------------------------#
# Plot density distribution of log-transformed FPKM values across samples
FPKM_Dist <- ggplot(fpkm_long, aes(x = log10(FPKM + 1), color = Sample_Name)) + 
  geom_density() +
  theme_classic() +
  labs(title = "Distribution of FPKM Values For Samples", x = "log10(FPKM + 1)", y = "Density")

# Display the plot
FPKM_Dist 

# Save the density plot
ggsave(paste0(fig_dir,"FPKM_Dist_plot.png"), plot = FPKM_Dist, width = 10, height = 6, dpi = 300)

# Create a violin plot of log-transformed FPKM values by sample
FPKM_Violin <- ggplot(fpkm_long, aes(x = Sample_Name, y = log10(FPKM + 1), fill = Sample_Name)) +
  geom_violin(trim = FALSE) +
  theme_classic() +
  labs(title = "FPKM Value Distribution (Violin)", x = "Sample", y = "log10(FPKM + 1)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# Save the violin plot
ggsave(paste0(fig_dir,"FPKM_Dist_ViolinPlot.png"), plot = FPKM_Violin, width = 10, height = 6, dpi = 300)

#----------------------------------#
# Heatmap of Top 60 Variable Genes
#----------------------------------#

# Add a log-transformed FPKM column
fpkm_long <- fpkm_long %>%
  mutate(logFPKM = log10(FPKM + 1))

# Pivot the long format back to wide format with genes as rows and samples as columns
fpkm_wide <- fpkm_long %>%
  select(GeneName, Sample_Name, logFPKM) %>%
  pivot_wider(
    names_from = Sample_Name,
    values_from = logFPKM,
    values_fn = mean
  )

# Drop rows with NA values and ensure unique gene names
fpkm_wide <- fpkm_wide %>%
  drop_na() %>%
  filter(!duplicated(GeneName))

# Convert wide-format FPKM data to a matrix for downstream analysis
fpkm_mat <- fpkm_wide %>%
  column_to_rownames("GeneName") %>%
  as.matrix()

# Compute variance for each gene across samples
var_genes <- apply(fpkm_mat, 1, var, na.rm = TRUE)

# Select the top 60 most variable genes that are also present in the matrix
top_genes <- intersect(names(sort(var_genes, decreasing = TRUE)[1:60]), rownames(fpkm_mat))
top_mat <- fpkm_mat[top_genes, ]

# Create a custom green-black-red color palette for the heatmap
green_red_palette <- colorRampPalette(c("green", "black", "red"))(100)

# Generate a heatmap of the top 60 variable genes
FPKMheatmap <- pheatmap(top_mat,
                        scale = "row",
                        cluster_rows = TRUE,
                        cluster_cols = TRUE,
                        color = green_red_palette,
                        fontsize_row = 6,
                        main = "Top 60 Most Variable Genes")

# Save the heatmap plot
ggsave(paste0(fig_dir,"FPKM_heatmap.png"), plot = FPKMheatmap, width = 10, height = 6, dpi = 300)

