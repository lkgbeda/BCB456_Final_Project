# This code determines the genes that are differentially expressed between Days.
#-----------------------------------------
# Load required libraries
#-----------------------------------------
library(DESeq2)
library(readr)
library(plyr)
library(dplyr)
library(tidyverse)
library(matrixStats)
library(stringr)
library(ggrepel)
library(pheatmap)
library(ggplot2)
library(VennDiagram)
library(grid)

# Define the base directory where your project files are located
base_dir <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project"
# Define the directory containing HTSeq count files
count_dir <- paste0(base_dir, "/HTSeq_Count/Cleaned_Count_Updated/")
fig_dir = paste0(base_dir,"/Figures/")
deg_file_dir = paste0(base_dir,"/DEA_DESeq/")

#--------------------------------------------------------------
# Load and Prepare Sample Metadata and HTSeq Count Files
#--------------------------------------------------------------

# List all HTSeq count files in the count directory that match the pattern "SRR<digits>_counts.txt"
sample_files <- list.files(count_dir, pattern = "SRR\\d+_counts\\.txt", full.names = FALSE)

# Preview the first few and last few count file names
head(sample_files)
tail(sample_files)

# Extract sample names from file names using regular expressions (e.g., 'SRR15969055')
sample_names <- stringr::str_extract(basename(sample_files), "SRR\\d+")
head(sample_names)
tail(sample_names)

# Import metadata table from a CSV file
sample_table <- read.table(
  paste0(base_dir,"/metadata.csv"),   # construct the full path to the metadata file
  header = TRUE,                      # use the first row as column names
  sep = ",",                          # fields are comma-separated
  fill = TRUE,                        # fill rows with unequal length
  quote = "\"",                       # allow quoted strings
  comment.char = "",                  # don't treat any character as comment
  stringsAsFactors = FALSE            # keep strings as characters, not factors
)

# Add the count file names as a new column in the metadata
sample_table$file_name <- sample_files
head(sample_table)

# Rename columns 1 and 4 to "SRX" and "SRR" for clarity
colnames(sample_table)[c(1,4)] <- c("SRX", "SRR")
head(sample_table)

# Reorder columns of the sample table (custom order likely for analysis compatibility)
sample_table <- sample_table[,c(2,8,7,1,6,5,2,3,4,6)]

# Sort metadata by the 'Replicate' column to maintain order
sample_table <- sample_table[order(sample_table$Replicate),]
head(sample_table)
tail(sample_table)

# Display the distribution of samples across Days and Groups (e.g., for experimental design validation)
table(sample_table$Day, sample_table$Group)

# Check if all count files exist (uncomment to run this validation)
# all(file.exists(file.path(directory, sample_table$sampleName)))

# Convert 'Day' and 'Group' columns to factors for use in modeling and analysis
sample_table$Day <- factor(sample_table$Day)
sample_table$Group <- factor(sample_table$Group)

# Display the levels of the 'Day' factor (e.g., Day 0, Day 3, Day 4, etc.)
levels(sample_table$Day)

# Display the levels of the 'Group' factor (e.g., Control, Treated, etc.)
levels(sample_table$Group)

#--------------------------------------------------------------
# Differential Expression Analysis by Condition (Day-based)
#--------------------------------------------------------------

# Create a new 'Condition' factor based on 'Day'
sample_table$Condition <- factor(sample_table$Day) 
# Create a DESeq2 dataset from HTSeq-count files using the sample table and specified directory
# Design formula specifies that differences will be modeled based on 'Condition'
Day_dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                      directory = count_dir,
                                      design = ~ Condition)

# Run the differential expression analysis pipeline on the constructed DESeq2 object
Day_dds <- DESeq(Day_dds)

#---------------------------------------------------------
## Transform Data for PCA
#----------------------------------------------------------

#----------------------------------------------------------
# Use variance stabilizing transformation (VST) to normalize counts before PCA
#----------------------------------------------------------
Day_vsd <- vst(Day_dds, blind = TRUE)
Day_gene_variances <- rowVars(assay(Day_dds))

# Select the top 200 most variable genes
Day_top_variable_genes <- order(Day_gene_variances, decreasing = TRUE)[1:2000]

# Subset VST-transformed matrix
Day_vsd_mat <- assay(Day_vsd)  # Extract transformed values
Day_vsd_top2000 <- Day_vsd_mat[Day_top_variable_genes, ]  # Keep only top 200 variable genes

Day_pca_result <- prcomp(t(Day_vsd_top2000), scale. = TRUE)  # Transpose to have samples as rows

# Extract PCA data for visualization
Day_pca_data <- as.data.frame(Day_pca_result$x)
Day_pca_data$Day <- colData(Day_vsd)$Day  # Add sample metadata
Day_pca_data$Group <- colData(Day_vsd)$Group
Day_vsd_pca <- ggplot(Day_pca_data, aes(PC1, PC2, color = Day,shape = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(100 * (Day_pca_result$sdev[1]^2 / sum(Day_pca_result$sdev^2)), 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * (Day_pca_result$sdev[2]^2 / sum(Day_pca_result$sdev^2)), 1), "% variance")) +
  ggtitle("PCA of Top 2000 Most Variable Genes (VST Transformation)") +
  theme_classic()

Day_vsd_pca
ggsave(paste0(fig_dir,"Day_vsd_pca_plot.png"), plot = Day_vsd_pca, width = 10, height = 6, dpi = 300)

#----------------------------------------------------------
## Transform Data for PCA using rlog
#----------------------------------------------------------

#rlog (regularized log transformation) to normalize counts before PCA
Day_rld <- rlog(Day_dds, blind = TRUE)
# Use rlog transformation instead of VST
Day_rld <- rlog(Day_dds, blind = TRUE)
Day_rld_gene_variances <- rowVars(assay(Day_rld))

# Select the top 200 most variable genes
Day_rld_top_variable_genes <- order(Day_rld_gene_variances, decreasing = TRUE)[1:200]

# Subset rlog-transformed matrix
Day_rld_mat <- assay(Day_rld)  # Extract transformed values
Day_rld_top2000 <- Day_rld_mat[Day_rld_top_variable_genes, ]  # Keep only top 200 variable genes

# Perform PCA
Day_rld_pca_result <- prcomp(t(Day_rld_top2000), scale. = TRUE)  # Transpose to have samples as rows

# Extract PCA data for visualization
Day_rld_pca_data <- as.data.frame(Day_rld_pca_result$x)
Day_rld_pca_data$Day <- colData(Day_rld)$Day# Add sample metadata
Day_rld_pca_data$Group <- colData(Day_rld)$Group # Add sample metadata

# Plot PCA with ggplot2
Day_rldpca <- ggplot(Day_rld_pca_data, aes(PC1, PC2, color = Day, shape = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(100 * (Day_rld_pca_result$sdev[1]^2 / sum(Day_rld_pca_result$sdev^2)), 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * (Day_rld_pca_result$sdev[2]^2 / sum(Day_rld_pca_result$sdev^2)), 1), "% variance")) +
  ggtitle("PCA of Top 2000 Most Variable Genes (rlog Transformation)") +
  theme_classic()  # Use a clean theme
Day_rldpca
ggsave(paste0(fig_dir,"Day_rldpca_plot.png"), plot = Day_rldpca, width = 10, height = 6, dpi = 300)

#----------------------------------------------------------
# Estimate dispersion
#----------------------------------------------------------
plotDispEsts(Day_dds)
png(paste0(fig_dir,"Day_dispersion_plot.png"), width = 800, height = 600)
plotDispEsts(Day_dds)
dev.off()

#--------------------------------------------------------------------------------
# Differential Expression Analysis and Gene Annotation for Time Point Comparisons 

#This script performs pairwise differential expression analysis using DESeq2
# for samples collected at different days (Day 0, Day 3, Day 4). The results 
# are annotated with gene names from an external mapping file, and formatted 
# for easier downstream use.
#--------------------------------------------------------------------------------
# Get results for each  Day vs The other day 
# Get results for each Day vs another Day comparison
res_day0_day3 <- results(Day_dds, contrast = c("Condition", "0", "3"))  # Day 0 vs Day 3
res_day0_day4 <- results(Day_dds, contrast = c("Condition", "0", "4"))  # Day 0 vs Day 4
res_day3_day4 <- results(Day_dds, contrast = c("Condition", "3", "4"))  # Day 3 vs Day 4

# Convert DESeqResults objects to data frames for downstream processing
res_day0_day3_df <- as.data.frame(res_day0_day3)
res_day0_day4_df <- as.data.frame(res_day0_day4)
res_day3_day4_df <- as.data.frame(res_day3_day4)

# Load gene name mapping file (contains gene IDs and their external names)
gene_name_mapping <- read.delim(paste0(base_dir,"/new_geneinfo.tsv"), sep = "\t", stringsAsFactors = FALSE)

# Add GeneID column to each data frame from row names (if not already present)
res_day0_day3_df$GeneID <- rownames(res_day0_day3_df)
res_day0_day4_df$GeneID <- rownames(res_day0_day4_df)
res_day3_day4_df$GeneID <- rownames(res_day3_day4_df)

# Join external gene names by GeneID
res_day0_day3_df <- res_day0_day3_df %>% left_join(gene_name_mapping, by = "GeneID")
res_day0_day4_df <- res_day0_day4_df %>% left_join(gene_name_mapping, by = "GeneID")
res_day3_day4_df <- res_day3_day4_df %>% left_join(gene_name_mapping, by = "GeneID")

# Replace missing 'old' gene names with GeneID
na_idx_03 <- is.na(res_day0_day3_df$old)
res_day0_day3_df$old[na_idx_03] <- res_day0_day3_df$GeneID[na_idx_03]
na_idx_04 <- is.na(res_day0_day4_df$old)
res_day0_day4_df$old[na_idx_04] <- res_day0_day4_df$GeneID[na_idx_04]
na_idx_34 <- is.na(res_day3_day4_df$old)
res_day3_day4_df$old[na_idx_34] <- res_day3_day4_df$GeneID[na_idx_34]

# Replace missing 'Name' column values with GeneID
na_idx_03 <- is.na(res_day0_day3_df$Name)
res_day0_day3_df$Name[na_idx_03] <- res_day0_day3_df$GeneID[na_idx_03]
na_idx_04 <- is.na(res_day0_day4_df$Name)
res_day0_day4_df$Name[na_idx_04] <- res_day0_day4_df$GeneID[na_idx_04]
na_idx_34 <- is.na(res_day3_day4_df$Name)
res_day3_day4_df$Name[na_idx_34] <- res_day3_day4_df$GeneID[na_idx_34]

# Function to reorder columns: puts Name, GeneID, and old first, followed by original results columns
reorder_columns <- function(df) {
  cols <- colnames(df)
  new_order <- c("Name", "GeneID", "old", cols[1:3], setdiff(cols, c("Name", "GeneID", "old", cols[1:3])))
  df <- df[, new_order]
  return(df)
}

# Apply column reordering to all results data frames
res_day0_day3_df <- reorder_columns(res_day0_day3_df)
res_day0_day4_df <- reorder_columns(res_day0_day4_df)
res_day3_day4_df <- reorder_columns(res_day3_day4_df)

# Preview the final data frames
head(res_day0_day3_df)
head(res_day0_day4_df)
head(res_day3_day4_df)


# ===============================
# Identify and Export Differentially Expressed Genes (Up- and Down-Regulated)
# ===============================

# -----------------------------------------------------------------------
# Genes Upregulated Between Days (log2FC > 1 and adjusted p-value < 0.05)
# -----------------------------------------------------------------------

up_Day0_vs_Day3 <- filter(res_day0_day3_df, padj < 0.05, log2FoldChange > 1)
write.table(up_Day0_vs_Day3, file = paste0(deg_file_dir, "up_Day0_vs_Day3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
head(up_Day0_vs_Day3)
up_Day0_vs_Day3_top60 <- up_Day0_vs_Day3 %>% arrange(padj) %>% slice_head(n = 60)

up_Day0_vs_Day4 <- filter(res_day0_day4_df, padj < 0.05, log2FoldChange > 1)
write.table(up_Day0_vs_Day4, file = paste0(deg_file_dir, "up_Day0_vs_Day4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
up_Day0_vs_Day4_top60 <- up_Day0_vs_Day4 %>% arrange(padj) %>% slice_head(n = 60)

up_Day3_vs_Day4 <- filter(res_day3_day4_df, padj < 0.05, log2FoldChange > 1)
write.table(up_Day3_vs_Day4, file = paste0(deg_file_dir, "up_Day3_vs_Day4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
up_Day3_vs_Day4_top60 <- up_Day3_vs_Day4 %>% arrange(padj) %>% slice_head(n = 60)


# --------------------------------------------------------------------------
# Genes Downregulated Between Days (log2FC < -1 and adjusted p-value < 0.05)
# --------------------------------------------------------------------------

down_Day0_vs_Day3 <- filter(res_day0_day3_df, padj < 0.05, log2FoldChange < -1)
write.table(down_Day0_vs_Day3, file = paste0(deg_file_dir, "down_Day0_vs_Day3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
down_Day0_vs_Day3_top60 <- down_Day0_vs_Day3 %>% arrange(padj) %>% slice_head(n = 60)

down_Day0_vs_Day4 <- filter(res_day0_day4_df, padj < 0.05, log2FoldChange < -1)
write.table(down_Day0_vs_Day4, file = paste0(deg_file_dir, "down_Day0_vs_Day4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
down_Day0_vs_Day4_top60 <- down_Day0_vs_Day4 %>% arrange(padj) %>% slice_head(n = 60)

down_Day3_vs_Day4 <- filter(res_day3_day4_df, padj < 0.05, log2FoldChange < -1)
write.table(down_Day3_vs_Day4, file = paste0(deg_file_dir, "down_Day3_vs_Day4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
down_Day3_vs_Day4_top60 <- down_Day3_vs_Day4 %>% arrange(padj) %>% slice_head(n = 60)

# ========================================================================
# Normalized Count Extraction and Heatmap Visualization of Selected Genes
# ========================================================================

# Get normalized counts using variance stabilizing transformation (VST)
Day_vsd <- vst(Day_dds, blind = FALSE)
Day_norm_counts <- assay(Day_vsd)
head(Day_norm_counts)
tail(Day_norm_counts)

# Define green-black-red color palette for heatmaps
green_red_palette <- colorRampPalette(c("green", "black", "red"))(100)

# Function to plot a heatmap for a list of genes
plot_heatmap <- function(gene_list, title) {
  # Get gene IDs and corresponding external names
  gene_ids <- gene_list$GeneID
  names <- gene_list$Name
  
  # Keep only gene IDs that are present in normalized counts
  valid_indices <- gene_ids %in% rownames(Day_norm_counts)
  gene_ids <- gene_ids[valid_indices]
  names <- names[valid_indices]
  
  # Subset the normalized counts matrix
  selected_counts <- Day_norm_counts[gene_ids, , drop = FALSE]
  
  # Set external names as rownames
  rownames(selected_counts) <- names
  
  # Z-score normalization per gene
  scaled_counts <- t(scale(t(selected_counts)))
  
  # Clip extreme values to improve color scale contrast
  scaled_counts[scaled_counts > 2] <- 2
  scaled_counts[scaled_counts < -2] <- -2
  
  # Optionally select top variable genes (e.g., top 50–60)
  top_genes <- rownames(scaled_counts)[1:min(60, nrow(scaled_counts))]
  plot_scaled_counts <- scaled_counts[top_genes, , drop = FALSE]
  
  # Plot the heatmap
  p <- pheatmap(plot_scaled_counts,
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                show_rownames = TRUE,
                show_colnames = TRUE,
                color = green_red_palette,
                breaks = seq(-2, 2, length.out = 101),
                main = title,
                fontsize_row = 6, 
                border_color = NA)
  
  print(p)
  return(p)
}


# ================================================================================
# Visualize Differential Expression of Selected Genes of Interest (from S3 and S4 Paper Data)
# ================================================================================

# -----------------------------------------------------
# Load Gene Lists from Supplementary Tables
# -----------------------------------------------------

# Read gene IDs from supplementary data files (S3 and S4)
df_geneIDs_S3 <- read.table(paste0(base_dir, "/all_geneIDs_S3.tsv"), sep = "\t", header = TRUE)
df_geneIDs_S4 <- read.table(paste0(base_dir, "/all_geneIDs_S4.tsv"), sep = "\t", header = TRUE)

# Check overlap between the two gene sets and combine unique gene entries
length(intersect(df_geneIDs_S3$old, df_geneIDs_S4$old))
geneIDs_S3_S4 <- unique(rbind(df_geneIDs_S3, df_geneIDs_S4)) 


# --------------------------------------------------------------
# Filter DEGs to Retain Only Genes of Interest for Visualization
# --------------------------------------------------------------

# For Day0 vs Day3, keep full up/down lists (not filtered by S3/S4 for now)
up_Day0_vs_Day3_Selected <- up_Day0_vs_Day3 
down_Day0_vs_Day3_Selected <- down_Day0_vs_Day3

# Filter Day0 vs Day4 and Day3 vs Day4 up/downregulated DEGs to retain only genes from S3/S4
up_Day0_vs_Day4_Selected <- up_Day0_vs_Day4 %>% filter(up_Day0_vs_Day4[[3]] %in% geneIDs_S3_S4$old) 
down_Day0_vs_Day4_Selected <- down_Day0_vs_Day4 %>% filter(down_Day0_vs_Day4[[3]] %in% geneIDs_S3_S4$old) 

up_Day3_vs_Day4_Selected <- up_Day3_vs_Day4 %>% filter(up_Day3_vs_Day4[[3]] %in% geneIDs_S3_S4$old) 
down_Day3_vs_Day4_Selected <- down_Day3_vs_Day4 %>% filter(down_Day3_vs_Day4[[3]] %in% geneIDs_S3_S4$old) 


# -----------------------------------------------------------
# Generate and Display Heatmaps for Each Timepoint Comparison
# -----------------------------------------------------------

# Day 0 vs Day 3
print(plot_heatmap(up_Day0_vs_Day3_Selected , "Genes Upregulated in Day0 vs Day3"))
print(plot_heatmap(down_Day0_vs_Day3_Selected, "Genes Downregulated in Day0 vs Day3"))

# Day 0 vs Day 4
print(plot_heatmap(up_Day0_vs_Day4 , "Genes Upregulated in Day0 vs Day4"))
print(plot_heatmap(down_Day0_vs_Day4, "Genes Downregulated in Day0 vs Day4"))

# Day 3 vs Day 4
print(plot_heatmap(up_Day3_vs_Day4, "Genes Upregulated in Day3 vs Day4"))
print(plot_heatmap(down_Day3_vs_Day4, "Genes Downregulated in Day3 vs Day4"))


# ----------------------------------------------------------------
# Save Heatmaps to figure directory as PNG Files
# ----------------------------------------------------------------

# Day 0 vs Day 3
ggsave(filename = paste0(fig_dir, "Heatmap_Day0_vs_Day3_upregulated.png"),
       plot = plot_heatmap(up_Day0_vs_Day3_Selected, "Genes Upregulated in Day0 vs Day3"),
       width = 6, height = 8)

ggsave(filename = paste0(fig_dir, "Heatmap_Day0_vs_Day3_downregulated.png"),
       plot = plot_heatmap(down_Day0_vs_Day3_Selected, "Genes Downregulated in Day0 vs Day3"),
       width = 6, height = 8)

# Day 0 vs Day 4
ggsave(filename = paste0(fig_dir, "Heatmap_Day0_vs_Day4_upregulated.png"),
       plot = plot_heatmap(up_Day0_vs_Day4, "Genes Upregulated in Day0 vs Day4"),
       width = 6, height = 8)

ggsave(filename = paste0(fig_dir, "Heatmap_Day0_vs_Day4_downregulated.png"),
       plot = plot_heatmap(down_Day0_vs_Day4, "Genes Downregulated in Day0 vs Day4"),
       width = 6, height = 8)

# Day 3 vs Day 4
ggsave(filename = paste0(fig_dir, "Heatmap_Day3_vs_Day4_upregulated.png"),
       plot = plot_heatmap(up_Day3_vs_Day4, "Genes Upregulated in Day3 vs Day4"),
       width = 6, height = 8)

ggsave(filename = paste0(fig_dir, "Heatmap_Day3_vs_Day4_downregulated.png"),
       plot = plot_heatmap(down_Day3_vs_Day4, "Genes Downregulated in Day3 vs Day4"),
       width = 6, height = 8)


# ===============================
# Summary Bar Plot of Gene Expression Regulation Across Timepoint Comparisons
# ===============================

# --------------------------------
# Count Total Genes per Comparison (from full DESeq2 results)
# --------------------------------
total_Day0_vs_Day3 <- nrow(res_day0_day3_df)
total_Day0_vs_Day4 <- nrow(res_day0_day4_df)
total_Day3_vs_Day4 <- nrow(res_day3_day4_df)

# --------------------------------
# Count Upregulated, Downregulated, and Nonregulated Genes
# --------------------------------
up_counts <- c(nrow(up_Day0_vs_Day3), nrow(up_Day0_vs_Day4), nrow(up_Day3_vs_Day4))
down_counts <- c(nrow(down_Day0_vs_Day3), nrow(down_Day0_vs_Day4), nrow(down_Day3_vs_Day4))

# Calculate non-regulated genes (not significantly up or down)
non_counts <- c(
  total_Day0_vs_Day3 - up_counts[1] - down_counts[1],
  total_Day0_vs_Day4 - up_counts[2] - down_counts[2],
  total_Day3_vs_Day4 - up_counts[3] - down_counts[3]
)

# --------------------------------
# Construct Data Frame for Plotting
# --------------------------------
comparison <- c("Day0 vs Day3", "Day0 vs Day4", "Day3 vs Day4")

df_plot <- data.frame(
  Comparison = rep(comparison, each = 3),
  Regulation = rep(c("Upregulated", "Downregulated", "Nonregulated"), times = 3),
  Count = c(up_counts, down_counts, non_counts)
)

# Ensure desired order for regulation types in the legend
df_plot <- df_plot %>%
  mutate(Regulation = factor(Regulation, levels = c("Upregulated", "Downregulated", "Nonregulated")))

# --------------------------------
# Create Bar Plot with Regulation Status Per Timepoint Comparison
# --------------------------------
Day_Regulation_bar <- ggplot(df_plot, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, size = 3.5) +
  labs(
    title = "Gene Expression Summary Across Days",
    y = "Number of Differentially Expressed Genes",
    x = ""
  ) +
  scale_fill_manual(
    values = c("Upregulated" = "#B22222", 
               "Downregulated" = "#4682B4", 
               "Nonregulated" = "#9370DB")
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# --------------------------------
# Save the Bar Plot as PNG
# --------------------------------
ggsave(filename = paste0(fig_dir, "Day_gene_expression_summary.png"), 
       plot = Day_Regulation_bar, 
       width = 10, height = 6)
