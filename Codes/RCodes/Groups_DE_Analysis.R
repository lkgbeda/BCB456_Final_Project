# This code determines the genes that are differentially expressed between CD and ND.
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
library(ggVennDiagram)
library(grid)


# ===============================================================
#  Define Project Directories where project files are located
# ===============================================================

# Base directory for the project
base_dir <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project"

# Directory containing cleaned HTSeq count files
count_dir <- paste0(base_dir, "/HTSeq_Count/Cleaned_Count_Updated/")

# Directory to save figures
fig_dir <- paste0(base_dir, "/Figures/")

# Directory for differential expression analysis results
deg_file_dir <- paste0(base_dir, "/DEA_DESeq/")

# (Commented out) Hardcoded alternatives to the above dynamically constructed paths
# directory <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/HTSeq_Count/Cleaned_Count_Updated/"
# fig_dir <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Figures/"
# deg_file_dir <- "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/DEA_DESeq/"

# ======================================================
#  Identify and Extract Sample Files
# ======================================================

# List all HTSeq count files with names matching pattern 'SRR######_counts.txt'
sample_files <- list.files(count_dir, pattern = "SRR\\d+_counts\\.txt", full.names = FALSE)

# Preview the first and last few filenames
head(sample_files)
tail(sample_files)

# Extract sample names (e.g., "SRR15969055") using regex
sample_names <- stringr::str_extract(basename(sample_files), "SRR\\d+")

# Preview extracted sample names
head(sample_names)
tail(sample_names)

# ====================================================
# Load and Prepare Sample Table
# ====================================================

# Read metadata CSV containing sample information
sample_table <- read.table(
  paste0(base_dir, "/metadata.csv"),  # Full path to metadata file
  header = TRUE,                      # Use the first row as header
  sep = ",",                          # CSV file (comma-separated)
  fill = TRUE,                        # Fill missing values in rows
  quote = "\"",                       # Allow strings in quotes
  comment.char = "",                  # Disable comment character
  stringsAsFactors = FALSE            # Keep character columns as strings
)

# Add file names to the metadata
sample_table$file_name <- sample_files

# Preview the updated metadata
head(sample_table)

# Rename specific columns for consistency
colnames(sample_table)[c(1,4)] <- c("SRX", "SRR")

# Preview the renamed columns
head(sample_table)

# Reorder columns to desired format
sample_table <- sample_table[, c(2,8,7,1,6,5,2,3,4,6)]

# Sort metadata by Replicate number
sample_table <- sample_table[order(sample_table$Replicate), ]

# Preview sorted metadata
head(sample_table)
tail(sample_table)

# Display sample distribution across Day and Group
table(sample_table$Day, sample_table$Group)

# =======================================================
#  Format Factor Variables for DESeq2
# =======================================================

# Convert Day and Group columns to factors (important for DESeq2 analysis)
sample_table$Day <- factor(sample_table$Day)
sample_table$Group <- factor(sample_table$Group)

# Check the factor levels for verification
levels(sample_table$Day)
levels(sample_table$Group)


# ==========================================================
# 5. Create DESeq2 Dataset and Run DEA
# ==========================================================

# Create a new factor variable 'Condition' from the 'Group' column
# This is used as the experimental design factor in DESeq2
sample_table$Condition <- factor(sample_table$Group)  

# Construct DESeq2 dataset object from HTSeq count files and sample metadata
dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sample_table,    # Metadata table with condition information
  directory = count_dir,         # Directory containing HTSeq count files
  design = ~ Condition           # Experimental design formula
)

#-----------------------------------------------------------------
# Run the DESeq2 differential expression analysis pipeline
#-----------------------------------------------------------------
dds <- DESeq(dds)



# ==========================================
# Transform Data and Perform PCA Analysis
# ==========================================

#------------------------------------------------------------------------
# Use variance stabilizing transformation (VST) to normalize count data
#------------------------------------------------------------------------
# This transformation helps stabilize variance across the range of counts
vsd <- vst(dds, blind = TRUE)

# Compute variance for each gene across all samples
vsd_gene_variances <- rowVars(assay(vsd))

# Select the top 2000 most variable genes based on variance
top_variable_genes <- order(vsd_gene_variances, decreasing = TRUE)[1:2000]

# Subset the VST-transformed expression matrix to only the top variable genes
vsd_mat <- assay(vsd)  # Extract transformed expression matrix
vsd_top2000 <- vsd_mat[top_variable_genes, ]  # Keep only top 2000 variable genes

# Perform Principal Component Analysis (PCA)
# Transpose the matrix so that samples are rows and genes are columns
vsd_pca_result <- prcomp(t(vsd_top2000), scale. = TRUE)

# Extract PCA scores and convert to data frame for plotting
vsd_pca_data <- as.data.frame(vsd_pca_result$x)
vsd_pca_data$Group <- colData(vsd)$Group  # Add group information from sample metadata

# Create PCA scatter plot of PC1 vs PC2
vsd_pca <- ggplot(vsd_pca_data, aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(100 * (vsd_pca_result$sdev[1]^2 / sum(vsd_pca_result$sdev^2)), 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * (vsd_pca_result$sdev[2]^2 / sum(vsd_pca_result$sdev^2)), 1), "% variance")) +
  ggtitle("PCA of Top 2000 Most Variable Genes (VST Transformation)") +
  theme_classic()

# Display PCA plot
vsd_pca

# Save PCA plot to file
ggsave(paste0(fig_dir, "vsd_pca_plot.png"), plot = vsd_pca, width = 10, height = 6, dpi = 300)

# Compute proportion of variance explained by each principal component
vsd_variance_explained <- (vsd_pca_result$sdev^2) / sum(vsd_pca_result$sdev^2)

# Create data frame for elbow plot
elbow_data <- data.frame(
  PC = seq_along(vsd_variance_explained),
  Variance = vsd_variance_explained
)

# Generate elbow plot to visualize how much variance each PC explains
ggplot(elbow_data, aes(x = PC, y = Variance)) +
  geom_point(size = 3, color = "blue") +
  geom_line(group = 1, color = "blue") +
  xlab("Principal Component") +
  ylab("Proportion of Variance Explained") +
  ggtitle("Elbow Plot for PCA (VST Transformation)") +
  theme_minimal()

#-------------------------------------------------------------------------
## Transform Data for PCA using rlog
#-------------------------------------------------------------------------
#rlog (regularized log transformation) to normalize counts before PCA
rld <- rlog(dds, blind = TRUE)
# Use rlog transformation instead of VST
rld <- rlog(dds, blind = TRUE)
gene_variances <- rowVars(assay(rld))

# Select the top 2000 most variable genes
top_variable_genes <- order(gene_variances, decreasing = TRUE)[1:200]

# Subset rlog-transformed matrix
rld_mat <- assay(rld)  # Extract transformed values
rld_top2000 <- rld_mat[top_variable_genes, ]  # Keep only top 2000 variable genes

# Perform PCA
rld_pca_result <- prcomp(t(rld_top2000), scale. = TRUE)  # Transpose to have samples as rows

# Extract PCA data for visualization
rld_pca_data <- as.data.frame(rld_pca_result$x)
rld_pca_data$Group <- colData(rld)$Group # Add sample metadata

# Plot PCA with ggplot2
rldpca <- ggplot(rld_pca_data, aes(PC1, PC2, color = Group, shape = Group)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", round(100 * (rld_pca_result$sdev[1]^2 / sum(rld_pca_result$sdev^2)), 1), "% variance")) +
  ylab(paste0("PC2: ", round(100 * (rld_pca_result$sdev[2]^2 / sum(rld_pca_result$sdev^2)), 1), "% variance")) +
  ggtitle("PCA of Top 2000 Most Variable Genes (rlog Transformation)") +
  theme_classic()  # Use a clean theme
# Display PCA plot
rldpca
# Save PCA plot to file
ggsave(paste0(fig_dir,"rldpca_plot.png"), plot = rldpca, width = 10, height = 6, dpi = 300)

# Compute variance explained
rld_variance_explained <- (rld_pca_result$sdev^2) / sum(rld_pca_result$sdev^2)

# Convert to a data frame for plotting
rld_elbow_data <- data.frame(
  PC = seq_along(rld_variance_explained),
  Variance = rld_variance_explained
)

# Plot the elbow plot
ggplot(rld_elbow_data, aes(x = PC, y = Variance)) +
  geom_point(size = 3, color = "blue") +
  geom_line(group = 1, color = "blue") +
  xlab("Principal Component") +
  ylab("Proportion of Variance Explained") +
  ggtitle("Elbow Plot for PCA (rlog Transformation)") +
  theme_minimal()


# -------------------------------------------------------------------------
#Estimate and Plot Dispersion
# -------------------------------------------------------------------------

# Plot dispersion estimates for DESeq2 model (diagnostic plot)
plotDispEsts(dds)

# Save the dispersion plot to file as PNG
png(paste0(fig_dir, "dispersion_plot.png"), width = 800, height = 600)
plotDispEsts(dds)
dev.off()

#==========================================================================
# Differential Gene Expression: CD vs ND at 3 Timepoints
### Using DESeq2 results and gene name annotation
#==========================================================================

#-------------------------------------------------------------------------
#Get DE Results from DESeq2 for Each group 
#-------------------------------------------------------------------------

# Extract differential expression results comparing tea (CD) vs tobacco (ND)
# at each of three time points (Day 0, Day 3, and Day 4)
res_day0 <- results(dds, contrast = c("Condition", "CD0", "ND0"))  # Day 0 comparison
res_day3 <- results(dds, contrast = c("Condition", "CD3", "ND3"))  # Day 3 comparison
res_day4 <- results(dds, contrast = c("Condition", "CD4", "ND4"))  # Day 4 comparison

# Convert the DESeqResults objects into standard data frames for easier manipulation
res_day0_df <- as.data.frame(res_day0)
res_day3_df <- as.data.frame(res_day3)
res_day4_df <- as.data.frame(res_day4)


#================================================================
# Annotate Results with External Gene Names 
#================================================================

# Load external gene name annotations
# The mapping file should contain at least two columns: GeneID and gene name(s)
gene_name_mapping <- read.delim(paste0(base_dir,"/new_geneinfo.tsv"), sep = "\t", stringsAsFactors = FALSE)

# Add a column called 'GeneID' to each DESeq2 result by copying row names
# This will be used as the key to join with the gene name mapping file
res_day0_df$GeneID <- rownames(res_day0_df)
res_day3_df$GeneID <- rownames(res_day3_df)
res_day4_df$GeneID <- rownames(res_day4_df)

# Merge the DE results with the gene name mapping using a left join on GeneID
# This adds external gene name columns like 'Name' and 'old' to the result
res_day0_df <- res_day0_df %>% left_join(gene_name_mapping, by = "GeneID")
res_day3_df <- res_day3_df %>% left_join(gene_name_mapping, by = "GeneID")
res_day4_df <- res_day4_df %>% left_join(gene_name_mapping, by = "GeneID")

# Count how many entries are missing external names ('Name' or 'old')
# These are typically NA because they could not be matched during the join
sum(is.na(res_day0_df$name))
sum(is.na(res_day3_df$name))
sum(is.na(res_day4_df$name))

sum(is.na(res_day0_df$old))
sum(is.na(res_day3_df$old))
sum(is.na(res_day4_df$old))

sum(is.na(res_day0_df$GeneID))  # Should be 0 if rownames were added properly
sum(is.na(res_day3_df$GeneID))
sum(is.na(res_day4_df$GeneID))

# Fill in missing 'old' name entries with the GeneID (useful as fallback identifiers)
na_idx_day0 <- is.na(res_day0_df$old)
res_day0_df$old[na_idx_day0] <- res_day0_df$GeneID[na_idx_day0]

na_idx_day3 <- is.na(res_day3_df$old)
res_day3_df$old[na_idx_day3] <- res_day3_df$GeneID[na_idx_day3]

na_idx_day4 <- is.na(res_day4_df$old)
res_day4_df$old[na_idx_day4] <- res_day4_df$GeneID[na_idx_day4]

# Similarly, fill in missing 'Name' fields with GeneID
# This ensures all rows have a usable name field for downstream filtering or plotting
na_idx_day0 <- is.na(res_day0_df$Name)
res_day0_df$Name[na_idx_day0] <- res_day0_df$GeneID[na_idx_day0]

na_idx_day3 <- is.na(res_day3_df$Name)
res_day3_df$Name[na_idx_day3] <- res_day3_df$GeneID[na_idx_day3]

na_idx_day4 <- is.na(res_day4_df$Name)
res_day4_df$Name[na_idx_day4] <- res_day4_df$GeneID[na_idx_day4]

# Function to reorder columns so that identifier and name fields appear first
# This improves the readability and consistency of output files
reorder_columns <- function(df) {
  cols <- colnames(df)
  new_order <- c("Name", "GeneID", "old", cols[1:3], setdiff(cols, c("Name", "GeneID", "old", cols[1:3])))
  df <- df[, new_order]
  return(df)
}

# Apply column reordering to each day's result
res_day0_df <- reorder_columns(res_day0_df)
res_day3_df <- reorder_columns(res_day3_df)
res_day4_df <- reorder_columns(res_day4_df)

# Preview the top few rows of each annotated DE result data frame
head(res_day0_df)
head(res_day3_df)
head(res_day4_df)


#================================================================================
# Identify and Export Upregulated Genes in CD 
#================================================================================
#----------------------------------------------------------------
# Identify genes significantly upregulated in tea (CD) compared to tobacco (ND)
# Criteria: adjusted p-value < 0.05 and log2FoldChange > 1 (strong upregulation)
#----------------------------------------------------------------
up_CD0_vs_ND0 <- filter(res_day0_df, padj < 0.05, log2FoldChange > 1)
# Write the list of upregulated genes on Day 3 to a tab-separated file. 
# Write the list of upregulated genes on Day 0 to a tab-separated file. 
# The file will be saved to the directory specified by 'deg_file_dir' and named "up_CD0_vs_ND0.tsv".
write.table(up_CD0_vs_ND0, file = paste0(deg_file_dir, "up_CD0_vs_ND0.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
up_CD0_vs_ND0_top60 <- up_CD0_vs_ND0 %>% arrange(padj) %>% slice_head(n = 60)

up_CD3_vs_ND3 <- filter(res_day3_df, padj < 0.05, log2FoldChange > 1)
# Write the list of upregulated genes on Day 3 to a tab-separated file. 
# The file will be saved to the directory specified by 'deg_file_dir' and named "up_CD3_vs_ND3.tsv".
write.table(up_CD3_vs_ND3, file = paste0(deg_file_dir, "up_CD3_vs_ND3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Top 60 upregulated genes ranked by significance
up_CD3_vs_ND3_top60 <- up_CD3_vs_ND3 %>% arrange(padj) %>% slice_head(n = 60)

up_CD4_vs_ND4 <- filter(res_day4_df, padj < 0.05, log2FoldChange > 1)
# Write the list of upregulated genes on Day 4 to a tab-separated file. 
# The file will be saved to the directory specified by 'deg_file_dir' and named "up_CD4_vs_ND4.tsv".
write.table(up_CD4_vs_ND4, file = paste0(deg_file_dir, "up_CD4_vs_ND4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Top 60 upregulated genes ranked by significance
up_CD4_vs_ND4_top60 <- up_CD4_vs_ND4 %>% arrange(padj) %>% slice_head(n = 60)


#=====================================================================================
# Identify and Export Downregulated Genes in CD 
#(i.e., Genes Upregulated in ND)                
#=====================================================================================

#------------------------------------------------------------------------------------
# Identify genes significantly downregulated in tea (CD), i.e., upregulated in ND
# Criteria: adjusted p-value < 0.05 and log2FoldChange < -1
#-------------------------------------------------------------------------------------
down_CD0_vs_ND0 <- filter(res_day0_df, padj < 0.05, log2FoldChange < -1)
# Write the list of downregulated genes on Day 0 to a tab-separated file. 
# The file will be saved to the directory specified by 'deg_file_dir' and named "down_CD0_vs_ND0.tsv".
write.table(down_CD0_vs_ND0, file = paste0(deg_file_dir, "down_CD0_vs_ND0.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Top 60 downregulated genes ranked by significance
down_CD0_vs_ND0_top60 <- down_CD0_vs_ND0 %>% arrange(padj) %>% slice_head(n = 60)

down_CD3_vs_ND3 <- filter(res_day3_df, padj < 0.05, log2FoldChange < -1)
# Write the list of downregulated genes on Day 3 to a tab-separated file. 
# The file will be saved to the directory specified by 'deg_file_dir' and named "down_CD3_vs_ND3.tsv".
write.table(down_CD3_vs_ND3, file = paste0(deg_file_dir, "down_CD3_vs_ND3.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Top 60 downregulated genes ranked by significance
down_CD3_vs_ND3_top60 <- down_CD3_vs_ND3 %>% arrange(padj) %>% slice_head(n = 60)

down_CD4_vs_ND4 <- filter(res_day4_df, padj < 0.05, log2FoldChange < -1) 
# Write the list of downregulated genes on Day 4 to a tab-separated file. 
# The file will be saved to the directory specified by 'deg_file_dir' and named "down_CD4_vs_ND4.tsv".
write.table(down_CD4_vs_ND4, file = paste0(deg_file_dir, "down_CD4_vs_ND4.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
# Top 60 downregulated genes ranked by significance
down_CD4_vs_ND4_top60 <- down_CD4_vs_ND4 %>% arrange(padj) %>% slice_head(n = 60)


#====================================================================
# UPREGULATED GENES CATEGORIZATION AND VISUALIZATION
#====================================================================

#---------------------------------------------
# Extract row names (i.e., gene IDs) of upregulated DEGs for each timepoint
#---------------------------------------------
genes_day0_up <- rownames(up_CD0_vs_ND0)   # Day 0 comparison
genes_day3_up <- rownames(up_CD3_vs_ND3)   # Day 3 comparison
genes_day4_up <- rownames(up_CD4_vs_ND4)   # Day 4 comparison

#---------------------------------------------
# Count the number of upregulated genes for each day
#---------------------------------------------
n_day0_up <- length(genes_day0_up)
n_day3_up <- length(genes_day3_up)
n_day4_up <- length(genes_day4_up)

#---------------------------------------------
# Calculate overlaps and exclusive gene sets
#---------------------------------------------
n_all_three_up <- length(Reduce(intersect, list(genes_day0_up, genes_day3_up, genes_day4_up)))  # Genes upregulated in all three comparisons

# Genes exclusively upregulated in one timepoint
n_only_day0_up <- length(setdiff(genes_day0_up, union(genes_day3_up, genes_day4_up)))
n_only_day3_up <- length(setdiff(genes_day3_up, union(genes_day0_up, genes_day4_up)))
n_only_day4_up <- length(setdiff(genes_day4_up, union(genes_day0_up, genes_day3_up)))

#---------------------------------------------
# Determine all genes tested and those not upregulated
#---------------------------------------------
all_genes <- union(union(rownames(res_day0_df), rownames(res_day3_df)), rownames(res_day4_df))  # All unique genes
all_up_DEGs <- union(union(genes_day0_up, genes_day3_up), genes_day4_up)                        # All upregulated genes
n_non_up_DEGs <- length(setdiff(all_genes, all_up_DEGs))                                        # Genes not significantly upregulated

#---------------------------------------------
# Pairwise upregulation overlaps (excluding third timepoint)
#---------------------------------------------
day0_day3_up_only <- setdiff(intersect(genes_day0_up, genes_day3_up), genes_day4_up)
day0_day4_up_only <- setdiff(intersect(genes_day0_up, genes_day4_up), genes_day3_up)
day3_day4_up_only <- setdiff(intersect(genes_day3_up, genes_day4_up), genes_day0_up)

#====================================================================
# PREPARE DATA FOR BAR PLOT
#====================================================================
bar_data_up <- data.frame(
  Category = c(
    "CD0_vs_ND0", "CD3_vs_ND3", "CD4_vs_ND4",                         # All DEGs
    "Only_CD0_vs_ND0", "Only_CD3_vs_ND3", "Only_CD4_vs_ND4",         # Exclusives
    "CD0_vs_ND0_and_CD3_vs_ND3", "CD0_vs_ND0_and_CD4_vs_ND4",        # Pairwise overlaps
    "CD3_vs_ND3_and_CD4_vs_ND4", "Shared_all_three", "Non_DEGs"      # Triple overlap and non-DEGs
  ),
  Count = c(
    n_day0_up, n_day3_up, n_day4_up,
    n_only_day0_up, n_only_day3_up, n_only_day4_up,
    length(day0_day3_up_only), length(day0_day4_up_only), length(day3_day4_up_only),
    n_all_three_up, n_non_up_DEGs
  )
)

#---------------------------------------------
# Generate bar plot
#---------------------------------------------
bar_plot_up <- ggplot(bar_data_up, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +                                 # Create bars
  geom_text(aes(label = Count), vjust = -0.3, size = 4) +       # Add text labels
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +    # Rotate labels
  labs(
    title = "Gene Expression Summary for Upregulated Genes (CD vs ND)",
    y = "Number of Differentially Expressed Genes", x = ""
  )

# Display the plot
print(bar_plot_up)

# Save bar plot to file
ggsave(paste0(fig_dir, "upregulated_gene_Bar.png"), plot = bar_plot_up, width = 10, height = 6, dpi = 300)

#====================================================================
# VENN DIAGRAM DATA CALCULATION
#====================================================================

# Recalculate exclusive and overlapping gene sets for clarity
only_day0_up <- setdiff(genes_day0_up, union(genes_day3_up, genes_day4_up))
only_day3_up <- setdiff(genes_day3_up, union(genes_day0_up, genes_day4_up))
only_day4_up <- setdiff(genes_day4_up, union(genes_day0_up, genes_day3_up))
all_three_up <- Reduce(intersect, list(genes_day0_up, genes_day3_up, genes_day4_up))
none_up <- setdiff(all_genes, all_up_DEGs)  # Genes not upregulated in any comparison

#---------------------------------------------
# Create a named list of counts for use in Venn plotting
#---------------------------------------------
venn_counts_up <- list(
  only_CD0 = length(only_day0_up),
  only_CD3 = length(only_day3_up),
  only_CD4 = length(only_day4_up),
  CD0_CD3 = length(day0_day3_up_only),
  CD0_CD4 = length(day0_day4_up_only),
  CD3_CD4 = length(day3_day4_up_only),
  all_three = length(all_three_up),
  none = length(none_up)
)

# Print counts for verification
print(venn_counts_up)
#====================================================================
# BASIC VENN DIAGRAM (STATIC - VennDiagram PACKAGE)
#====================================================================

# Load required packages
library(VennDiagram)
library(grid)

#---------------------------
# Prepare Venn Diagram
#---------------------------
venn.plot_up <- draw.triple.venn(
  area1 = length(only_day0_up) + length(day0_day3_up_only) + length(day0_day4_up_only) + length(all_three_up),
  area2 = length(only_day3_up) + length(day0_day3_up_only) + length(day3_day4_up_only) + length(all_three_up),
  area3 = length(only_day4_up) + length(day0_day4_up_only) + length(day3_day4_up_only) + length(all_three_up),
  n12 = length(day0_day3_up_only) + length(all_three_up),
  n23 = length(day3_day4_up_only) + length(all_three_up),
  n13 = length(day0_day4_up_only) + length(all_three_up),
  n123 = length(all_three_up),
  category = c("CD0 vs ND0", "CD3 vs ND3", "CD4 vs ND4"),
  fill = c("red", "orange", "purple"),
  cex = 2,
  cat.cex = 1.5,
  cat.col = c("red", "orange", "purple")
)

#---------------------------
# Draw to Screen (Interactive)
#---------------------------
grid.newpage()
grid.draw(venn.plot_up)
grid.text(paste("None:", length(none_up)), x = unit(0.9, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 15))

#---------------------------
# Save to PNG File
#---------------------------
png(filename = file.path(fig_dir, "upregulated_venn_diagram.png"), width = 600, height = 700)
grid.newpage()
pushViewport(viewport(layout = grid.layout(10, 1)))

# Add title
pushViewport(viewport(layout.pos.row = 1))
grid.text("Upregulated Gene Overlap", 
          x = unit(0.5, "npc"), y = unit(0.5, "npc"),
          gp = gpar(fontsize = 18, fontface = "bold"))
popViewport()

# Add Venn diagram
pushViewport(viewport(layout.pos.row = 2:9))
grid.draw(venn.plot_up)
popViewport()

# Add text for genes not included in any category
pushViewport(viewport(layout.pos.row = 10))
grid.text(paste("None:", length(none_up)), 
          x = unit(0.9, "npc"), y = unit(0.5, "npc"),
          gp = gpar(fontsize = 15))
popViewport()

# Finish PNG
dev.off()




#====================================================================
# BASIC VENN DIAGRAM (STATIC - VennDiagram PACKAGE)
#====================================================================
# Create static Venn diagram for upregulated genes
venn.plot_up <- draw.triple.venn(
  area1 = length(only_day0_up) + length(day0_day3_up_only) + length(day0_day4_up_only) + length(all_three_up),
  area2 = length(only_day3_up) + length(day0_day3_up_only) + length(day3_day4_up_only) + length(all_three_up),
  area3 = length(only_day4_up) + length(day0_day4_up_only) + length(day3_day4_up_only) + length(all_three_up),
  n12 = length(day0_day3_up_only) + length(all_three_up),
  n23 = length(day3_day4_up_only) + length(all_three_up),
  n13 = length(day0_day4_up_only) + length(all_three_up),
  n123 = length(all_three_up),
  category = c("CD0 vs ND0", "CD3 vs ND3", "CD4 vs ND4"),
  fill = c("red", "orange", "purple"),
  cex = 2,
  cat.cex = 1.5,
  cat.col = c("red", "orange", "purple")
  
)

# Draw on screen
grid.draw(venn.plot_up)

# Annotate genes not included in any DEG category
grid.text(paste("None:", length(none_up)), x = unit(0.9, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 15))

# Save the Venn diagram to a PNG file
png(paste0(fig_dir, "upregulated_venn_diagram.png"), width = 600, height = 700)
grid.newpage()
pushViewport(viewport(layout = grid.layout(10, 1)))

# Add title
grid.text("Upregulated Gene Overlap", 
          x = unit(0.5, "npc"), y = unit(0.5, "npc"),
          gp = gpar(fontsize = 18, fontface = "bold"),
          vp = viewport(layout.pos.row = 1))

# Add Venn diagram
pushViewport(viewport(layout.pos.row = 2:9))
grid.draw(venn.plot_up)
popViewport()

# Add text for non-upregulated genes
grid.text(paste("None:", length(none_up)), 
          x = unit(0.9, "npc"), y = unit(0.5, "npc"),
          gp = gpar(fontsize = 15),
          vp = viewport(layout.pos.row = 10))

# Finalize PNG
dev.off()

#====================================================================
# ggVennDiagram (More aesthetic visualization)
#====================================================================

# Optional: Generate dummy gene names for clarity (for example only)
only_CD0 <- paste0("gene_cd0_", 1:length(only_day0_up))
only_CD3 <- paste0("gene_cd3_", 1:length(only_day3_up))
only_CD4 <- paste0("gene_cd4_", 1:length(only_day4_up))
CD0_CD3 <- paste0("gene_cd0_cd3_", 1:length(day0_day3_up_only))
CD0_CD4 <- paste0("gene_cd0_cd4_", 1:length(day0_day4_up_only))
CD3_CD4 <- paste0("gene_cd3_cd4_", 1:length(day3_day4_up_only))
all_three <- paste0("gene_all_", 1:length(all_three_up))

# Use actual gene lists for ggVennDiagram plotting
genes_up_real <- list(
  "CD0 vs ND0" = genes_day0_up,
  "CD3 vs ND3" = genes_day3_up,
  "CD4 vs ND4" = genes_day4_up
)

# Create ggVennDiagram plot
p_up <- ggVennDiagram(genes_up_real, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "red") +
  ggtitle("Upregulated Genes Overlap") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )

# Display Venn diagram
grid.newpage()
print(p_up, vp = viewport(layout = grid.layout(1, 1)))
grid.text(paste("None:", length(none_up)), x = unit(0.70, "npc"), y = unit(0.05, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold"))

# Save final version to PNG
png(paste0(fig_dir, "upregulated_genes_venn.png"), width = 800, height = 900, res = 120)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 1)))
print(p_up, vp = viewport(layout.pos.row = 1))
grid.text(paste("None:", length(none_up)), x = unit(0.70, "npc"), y = unit(0.05, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

#====================================================================
# DOWNREGULATED GENES CATEGORIZATION  AND VISUALIZATION 
#===================================================================
#====================================================================
# DOWNREGULATED GENES CATEGORIZATION AND VISUALIZATION
#====================================================================

#---------------------------------------------
# Extract row names (i.e., gene IDs) of downregulated DEGs for each timepoint
#---------------------------------------------
genes_day0_down <- rownames(down_CD0_vs_ND0)  # Day 0 comparison
genes_day3_down <- rownames(down_CD3_vs_ND3)  # Day 3 comparison
genes_day4_down <- rownames(down_CD4_vs_ND4)  # Day 4 comparison

#---------------------------------------------
# Count the number of downregulated genes for each day
#---------------------------------------------
n_day0_down <- length(genes_day0_down)
n_day3_down <- length(genes_day3_down)
n_day4_down <- length(genes_day4_down)

#---------------------------------------------
# Calculate overlaps and exclusive gene sets
#---------------------------------------------
n_all_three_down <- length(Reduce(intersect, list(genes_day0_down, genes_day3_down, genes_day4_down)))  # Genes downregulated in all three conditions

# Genes downregulated only in one specific comparison
n_only_day0_down <- length(setdiff(genes_day0_down, union(genes_day3_down, genes_day4_down)))
n_only_day3_down <- length(setdiff(genes_day3_down, union(genes_day0_down, genes_day4_down)))
n_only_day4_down <- length(setdiff(genes_day4_down, union(genes_day0_down, genes_day3_down)))

#---------------------------------------------
# Determine all genes tested and those not downregulated
#---------------------------------------------
all_genes <- union(union(rownames(res_day0_df), rownames(res_day3_df)), rownames(res_day4_df))  # All unique genes across timepoints
all_down_DEGs <- union(union(genes_day0_down, genes_day3_down), genes_day4_down)               # All downregulated DEGs
n_non_down_DEGs <- length(setdiff(all_genes, all_down_DEGs))                                   # Genes not significantly downregulated

#---------------------------------------------
# Pairwise downregulation overlaps (excluding the third timepoint)
#---------------------------------------------
day0_day3_down_only <- setdiff(intersect(genes_day0_down, genes_day3_down), genes_day4_down)
day0_day4_down_only <- setdiff(intersect(genes_day0_down, genes_day4_down), genes_day3_down)
day3_day4_down_only <- setdiff(intersect(genes_day3_down, genes_day4_down), genes_day0_down)

#====================================================================
# PREPARE DATA FOR BAR PLOT
#====================================================================
bar_data_down <- data.frame(
  Category = c(
    "CD0_vs_ND0", "CD3_vs_ND3", "CD4_vs_ND4",                       # All DEGs by group
    "Only_CD0_vs_ND0", "Only_CD3_vs_ND3", "Only_CD4_vs_ND4",       # Exclusively downregulated in one condition
    "CD0_vs_ND0_and_CD3_vs_ND3", "CD0_vs_ND0_and_CD4_vs_ND4",      # Pairwise overlaps
    "CD3_vs_ND3_and_CD4_vs_ND4", "Shared_all_three", "Non_DEGs"    # Triple overlap and non-DEGs
  ),
  Count = c(
    n_day0_down, n_day3_down, n_day4_down,
    n_only_day0_down, n_only_day3_down, n_only_day4_down,
    length(day0_day3_down_only), length(day0_day4_down_only), length(day3_day4_down_only),
    n_all_three_down, n_non_down_DEGs
  )
)

#---------------------------------------------
# Generate bar plot
#---------------------------------------------
bar_plot_down <- ggplot(bar_data_down, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +                               # Draw bars
  geom_text(aes(label = Count), vjust = -0.3, size = 4) +     # Add count labels above bars
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for clarity
  labs(
    title = "Gene Expression Summary for Downregulated Genes (CD vs ND)",
    y = "Number of Differentially Expressed Genes", x = ""
  )

# Display plot
print(bar_plot_down)

# Save bar plot to file
ggsave(paste0(fig_dir, "downregulated_gene_Bar.png"), plot = bar_plot_down, width = 10, height = 6, dpi = 300)

#====================================================================
# VENN DIAGRAM DATA CALCULATION
#====================================================================

# Recalculate exclusives and overlaps for clarity and reuse
only_day0_down <- setdiff(genes_day0_down, union(genes_day3_down, genes_day4_down))
only_day3_down <- setdiff(genes_day3_down, union(genes_day0_down, genes_day4_down))
only_day4_down <- setdiff(genes_day4_down, union(genes_day0_down, genes_day3_down))
all_three_down <- Reduce(intersect, list(genes_day0_down, genes_day3_down, genes_day4_down))
none_down <- setdiff(all_genes, all_down_DEGs)  # Genes not downregulated in any comparison

#---------------------------------------------
# Create a list of counts for Venn diagram input
#---------------------------------------------
venn_counts_down <- list(
  only_CD0 = length(only_day0_down),
  only_CD3 = length(only_day3_down),
  only_CD4 = length(only_day4_down),
  CD0_CD3 = length(day0_day3_down_only),
  CD0_CD4 = length(day0_day4_down_only),
  CD3_CD4 = length(day3_day4_down_only),
  all_three = length(all_three_down),
  none = length(none_down)
)

# Print out the Venn counts for sanity check
print(venn_counts_down)

#====================================================================
# BASIC VENN DIAGRAM (STATIC - VennDiagram PACKAGE)
#====================================================================
library(VennDiagram)
library(grid)

# Draw basic triple Venn diagram using computed overlaps
venn.plot_down <- draw.triple.venn(
  area1 = length(only_day0_down) + length(day0_day3_down_only) + length(day0_day4_down_only) + length(all_three_down),
  area2 = length(only_day3_down) + length(day0_day3_down_only) + length(day3_day4_down_only) + length(all_three_down),
  area3 = length(only_day4_down) + length(day0_day4_down_only) + length(day3_day4_down_only) + length(all_three_down),
  n12 = length(day0_day3_down_only) + length(all_three_down),
  n23 = length(day3_day4_down_only) + length(all_three_down),
  n13 = length(day0_day4_down_only) + length(all_three_down),
  n123 = length(all_three_down),
  category = c("CD0 vs ND0", "CD3 vs ND3", "CD4 vs ND4"),
  fill = c("red", "green", "purple"),
  cex = 2,
  cat.cex = 1.5,
  cat.col = c("red", "green", "purple")
)

# Draw Venn diagram to active graphical device
grid.draw(venn.plot_down)

# Add count of non-downregulated genes outside Venn
grid.text(paste("None:", length(none_down)), x = unit(0.9, "npc"), y = unit(0.1, "npc"), gp = gpar(fontsize = 15))

# Save Venn diagram as image with better formatting
png(paste0(fig_dir, "downregulated_venn_diagram.png"), width = 600, height = 700)
grid.newpage()
pushViewport(viewport(layout = grid.layout(10, 1)))

# Add title
grid.text("Downregulated Gene Overlap", 
          x = unit(0.5, "npc"), y = unit(0.5, "npc"),
          gp = gpar(fontsize = 18, fontface = "bold"),
          vp = viewport(layout.pos.row = 1))

# Plot Venn diagram in center
pushViewport(viewport(layout.pos.row = 2:9))
grid.draw(venn.plot_down)
popViewport()

# Add annotation for genes not in any set
grid.text(paste("None:", length(none_down)), 
          x = unit(0.9, "npc"), y = unit(0.5, "npc"),
          gp = gpar(fontsize = 15),
          vp = viewport(layout.pos.row = 10))

# Close graphics device
dev.off()

#====================================================================
# ggVennDiagram (Better visualization, custom gene names)
#====================================================================
# Create illustrative dummy gene names for plot clarity (optional but recommended)
only_CD0 <- paste0("gene_cd0_", 1:length(only_day0_down))
only_CD3 <- paste0("gene_cd3_", 1:length(only_day3_down))
only_CD4 <- paste0("gene_cd4_", 1:length(only_day4_down))
CD0_CD3 <- paste0("gene_cd0_cd3_", 1:length(day0_day3_down_only))
CD0_CD4 <- paste0("gene_cd0_cd4_", 1:length(day0_day4_down_only))
CD3_CD4 <- paste0("gene_cd3_cd4_", 1:length(day3_day4_down_only))
all_three <- paste0("gene_all_", 1:length(all_three_down))

# Create actual gene sets for ggVennDiagram plotting
genes_down_real <- list(
  "CD0 vs ND0" = genes_day0_down,
  "CD3 vs ND3" = genes_day3_down,
  "CD4 vs ND4" = genes_day4_down
)

# Generate Venn diagram using ggVennDiagram
p_down <- ggVennDiagram(genes_down_real, label_alpha = 0, label = "count") +
  scale_fill_gradient(low = "white", high = "green") +
  ggtitle("Downregulated Genes Overlap") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none"
  )

# Display with grid
grid.newpage()
print(p_down, vp = viewport(layout = grid.layout(1, 1)))
grid.text(paste("None:", length(none_down)), x = unit(0.70, "npc"), y = unit(0.05, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold"))

# Save plot as PNG
png(paste0(fig_dir, "downregulated_genes_venn.png"), width = 800, height = 900, res = 120)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 1)))
print(p_down, vp = viewport(layout.pos.row = 1))
grid.text(paste("None:", length(none_down)), x = unit(0.70, "npc"), y = unit(0.05, "npc"),
          gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()


#======================================================================
#Count upregulated and downregulated genes
#======================================================================
# Upregulated counts
n_up_day0 <- nrow(up_CD0_vs_ND0)
n_up_day3 <- nrow(up_CD3_vs_ND3)
n_up_day4 <- nrow(up_CD4_vs_ND4)

# Downregulated counts (already defined in your code)
n_down_day0 <- nrow(down_CD0_vs_ND0)
n_down_day3 <- nrow(down_CD3_vs_ND3)
n_down_day4 <- nrow(down_CD4_vs_ND4)


# Total genes tested in each comparison
total_genes_day0 <- nrow(res_day0_df)
total_genes_day3 <- nrow(res_day3_df)
total_genes_day4 <- nrow(res_day4_df)

# Non-regulated = total - upregulated - downregulated
n_non_day0 <- total_genes_day0 - n_up_day0 - n_down_day0
n_non_day3 <- total_genes_day3 - n_up_day3 - n_down_day3
n_non_day4 <- total_genes_day4 - n_up_day4 - n_down_day4


# Create updated data frame with all three categories
bar_data_full <- data.frame(
  Comparison = rep(c("CD0 vs ND0", "CD3 vs ND3", "CD4 vs ND4"), each = 3),
  Regulation = rep(c("Upregulated", "Downregulated", "Nonregulated"), times = 3),
  Count = c(n_up_day0, n_down_day0, n_non_day0,
            n_up_day3, n_down_day3, n_non_day3,
            n_up_day4, n_down_day4, n_non_day4)
)

# Set desired order for Regulation
bar_data_full$Regulation <- factor(bar_data_full$Regulation,
                                   levels = c("Upregulated", "Downregulated", "Nonregulated"))
#-----------------------------------------------------------------------------
# Bar Plot to compare Regulations 
#-----------------------------------------------------------------------------
reg_bar_plot <- ggplot(bar_data_full, aes(x = Comparison, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Count), 
            vjust = -0.3, position = position_dodge(width = 0.8), size = 3) +
  theme_classic() +
  scale_fill_manual(values = c(
    "Upregulated" = "firebrick", 
    "Downregulated" = "steelblue", 
    "Nonregulated" = "#9370DB"
  )) +
  labs(
    title = "Gene Expression Summary Across CD vs ND",
    y = "Number of Defferencially Expressed Genes",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
  )

reg_bar_plot 
ggsave(paste0(fig_dir,"regulateion_Bar.png"), plot = reg_bar_plot, width = 10, height = 6, dpi = 300)


# Preview a few rows
# Top 10 upregulated and downregulated genes at Day 0
head(up_CD0_vs_ND0, 10)     # Up in CD0 vs ND0 (tea upregulated)
head(down_CD0_vs_ND0, 10)   # Down in CD0 vs ND0 (tobacco upregulated)

# Top 10 upregulated and downregulated genes at Day 3
head(up_CD3_vs_ND3, 10)     # Up in CD3 vs ND3
head(down_CD3_vs_ND3, 10)   # Down in CD3 vs ND3

# Top 10 upregulated and downregulated genes at Day 4
head(up_CD4_vs_ND4, 10)     # Up in CD4 vs ND4
head(down_CD4_vs_ND4, 10)   # Down in CD4 vs ND4

#=============================================================================
#             Heatmap Visualization of DEG Across Conditions            
#    (CD vs ND at Day 0, 3, and 4 — Upregulated and Downregulated)      
#=============================================================================

# Get normalized counts
vsd <- vst(dds, blind = FALSE)
norm_counts <- assay(vsd)
head(norm_counts)
tail(norm_counts)
# Define green-black-red color palette
green_red_palette <- colorRampPalette(c("green", "black", "red"))(100)


plot_heatmap <- function(gene_list, title) {
  # Get gene IDs and corresponding external names
  gene_ids <- gene_list$GeneID
  names <- gene_list$Name
  
  # Keep only gene IDs that are present in norm_counts
  valid_indices <- gene_ids %in% rownames(norm_counts)
  gene_ids <- gene_ids[valid_indices]
  names <- names[valid_indices]
  
  # Subset the normalized counts matrix
  selected_counts <- norm_counts[gene_ids, , drop = FALSE]
  
  # Set external names as rownames
  rownames(selected_counts) <- names
  
  # Get sample condition labels
  sample_conditions <- colData(vsd)$Condition
  condition_levels <- unique(sample_conditions)
  
  # Create a matrix of mean expression values by condition
  averaged_counts <- sapply(condition_levels, function(cond) {
    rowMeans(selected_counts[, sample_conditions == cond, drop = FALSE])
  })
  
  # Set row and column names
  #rownames(averaged_counts) <- names
  colnames(averaged_counts) <- condition_levels
  
  # Reorder columns manually
  custom_order <- c("ND0", "CD0", "ND3", "CD3", "ND4", "CD4")
  averaged_counts <- averaged_counts[, custom_order, drop = FALSE]
  
  # Z-score normalization per gene
  scaled_counts <- t(scale(t(averaged_counts)))
  
  # Clip extreme values
  scaled_counts[scaled_counts > 2] <- 2
  scaled_counts[scaled_counts < -2] <- -2
  
  # (Optional) Select top 50 most variable genes
  top_genes <- rownames(scaled_counts)[1:min(60, nrow(scaled_counts))]
  plot_scaled_counts <- scaled_counts[top_genes, , drop = FALSE]
  
  # Plot the heatmap
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
  
  print(p)   # <---- ADD this!
  return(p)  # (optional) if you want to capture the heatmap object
  
}

#Selected genes for Deg visulization based on the gene of interest from the paper 
df_geneIDs_S3 <-  read.table(paste0(base_dir,"/all_geneIDs_S3.tsv"), sep = "\t", header = TRUE)

df_geneIDs_S4 <- read.table(paste0(base_dir,"/all_geneIDs_S4.tsv"), sep = "\t", header = TRUE)

length(intersect(df_geneIDs_S3$old, df_geneIDs_S4$old))
geneIDs_S3_S4 <- unique(rbind(df_geneIDs_S3, df_geneIDs_S4)) 


up_CD0_vs_ND0_Selected <- up_CD0_vs_ND0 
down_CD0_vs_ND0_Selected <- down_CD0_vs_ND0
up_CD3_vs_ND3_Selected <- up_CD3_vs_ND3 %>% filter(up_CD3_vs_ND3[[3]] %in% geneIDs_S3_S4$old) 
down_CD3_vs_ND3_Selected <- down_CD3_vs_ND3 %>% filter(down_CD3_vs_ND3[[3]] %in% geneIDs_S3_S4$old) 
up_CD4_vs_ND4_Selected <- up_CD4_vs_ND4 %>% filter(up_CD4_vs_ND4[[3]] %in% geneIDs_S3_S4$old) 
down_CD4_vs_ND4_Selected <- down_CD4_vs_ND4 %>% filter(down_CD4_vs_ND4[[3]] %in% geneIDs_S3_S4$old) 


# Plot top 60 up- and downregulated genes for each timepoint
# Day 0
print(plot_heatmap(up_CD0_vs_ND0_Selected , "Upregulated in CD0 vs ND0"))
print(plot_heatmap(down_CD0_vs_ND0_Selected, "Downregulated in CD0 vs ND0"))

# Day 3
print(plot_heatmap(down_CD3_vs_ND3_Selected, "Upregulated in CD3 vs ND3"))
print(plot_heatmap(down_CD3_vs_ND3_Selected, "Downregulated in CD3 vs ND3"))

# Day 4
print(plot_heatmap(up_CD4_vs_ND4_Selected, "Upregulated in CD4 vs ND4"))
print(plot_heatmap(down_CD4_vs_ND4_Selected, "Downregulated in CD4 vs ND4"))

# Upregulated in CD0 vs ND0
ggsave(filename = paste0(fig_dir, "Heatmap_CD0_vs_ND0_upregulated.png"),
       plot_heatmap(up_CD0_vs_ND0_Selected, "Genes Upregulated in CD0 vs ND0"),
       width = 4, height = 6)

# Downregulated in CD0 vs ND0
ggsave(filename = paste0(fig_dir, "Heatmap_CD0_vs_ND0_downregulated.png"),
       plot_heatmap(down_CD0_vs_ND0_Selected, "Genes Downregulated in CD0 vs ND0"),
       width = 4, height = 6)

# Upregulated in CD3 vs ND3
ggsave(filename = paste0(fig_dir, "Heatmap_CD3_vs_ND3_upregulated.png"),
       plot_heatmap(up_CD3_vs_ND3_Selected, "Genes Upregulated in CD3 vs ND3"),
       width = 4, height = 6)

# Downregulated in CD3 vs ND3
ggsave(filename = paste0(fig_dir, "Heatmap_CD3_vs_ND3_downregulated.png"),
       plot_heatmap(down_CD3_vs_ND3_Selected, "Genes Downregulated in CD3 vs ND3"),
       width = 4, height = 6)

# Upregulated in CD4 vs ND4
ggsave(filename = paste0(fig_dir, "Heatmap_CD4_vs_ND4_upregulated.png"),
       plot_heatmap(up_CD4_vs_ND4_Selected, "Genes Upregulated in CD4 vs ND4"),
       width = 4, height = 6)

# Downregulated in CD4 vs ND4
ggsave(filename = paste0(fig_dir, "Heatmap_CD4_vs_ND4_downregulated.png"),
       plot_heatmap(down_CD3_vs_ND3_Selected, " Genes Downregulated in CD4 vs ND4"),
       width = 4, height = 6)



# Day 0
print(plot_heatmap(up_CD0_vs_ND0_top60, "Upregulated in CD0 vs ND0"))
print(plot_heatmap(down_CD0_vs_ND0_top60, "Downregulated in CD0 vs ND0"))

# Day 3
print(plot_heatmap(up_CD3_vs_ND3, "Upregulated in CD3 vs ND3"))
print(plot_heatmap(down_CD3_vs_ND3_top60, "Downregulated in CD3 vs ND3"))

# Day 4
print(plot_heatmap(up_CD4_vs_ND4_top60, "Top 60 Upregulated in CD4 vs ND4"))
print(plot_heatmap(down_CD4_vs_ND4_top60, "Top 60 Downregulated in CD4 vs ND4"))

# Optional: Plot full sets
print(plot_heatmap(up_CD0_vs_ND0, "All Upregulated in CD0 vs ND0"))
print(plot_heatmap(down_CD0_vs_ND0, "All Downregulated in CD0 vs ND0"))

print(plot_heatmap(up_CD3_vs_ND3, "All Upregulated in CD3 vs ND3"))
print(plot_heatmap(down_CD3_vs_ND3, "All Downregulated in CD3 vs ND3"))

print(plot_heatmap(up_CD4_vs_ND4, "All Upregulated in CD4 vs ND4"))
print(plot_heatmap(down_CD4_vs_ND4, "All Downregulated in CD4 vs ND4"))


# Save top 60 heatmaps

fig_dir = "/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Figures/"
# Upregulated in CD0 vs ND0
ggsave(filename = paste0(fig_dir, "heatmap_CD0_vs_ND0_upregulated.png"),
       plot_heatmap(up_CD0_vs_ND0_top60, "Top 60 Upregulated in CD0 vs ND0"),
       width = 8, height = 6)

# Downregulated in CD0 vs ND0
ggsave(filename = paste0(fig_dir, "heatmap_CD0_vs_ND0_downregulated.png"),
       plot_heatmap(down_CD0_vs_ND0_top60, "Top 60 Downregulated in CD0 vs ND0"),
       width = 8, height = 6)

# Upregulated in CD3 vs ND3
ggsave(filename = paste0(fig_dir, "heatmap_CD3_vs_ND3_upregulated.png"),
       plot_heatmap(up_CD3_vs_ND3_top60, "Top 60 Upregulated in CD3 vs ND3"),
       width = 8, height = 6)

# Downregulated in CD3 vs ND3
ggsave(filename = paste0(fig_dir, "heatmap_CD3_vs_ND3_downregulated.png"),
       plot_heatmap(down_CD3_vs_ND3_top60, "Top 60 Downregulated in CD3 vs ND3"),
       width = 8, height = 6)

# Upregulated in CD4 vs ND4
ggsave(filename = paste0(fig_dir, "heatmap_CD4_vs_ND4_upregulated.png"),
       plot_heatmap(up_CD4_vs_ND4_top60, "Top 60 Upregulated in CD4 vs ND4"),
       width = 8, height = 6)

# Downregulated in CD4 vs ND4
ggsave(filename = paste0(fig_dir, "heatmap_CD4_vs_ND4_downregulated.png"),
       plot_heatmap(down_CD4_vs_ND4_top60, "Top 60 Downregulated in CD4 vs ND4"),
       width = 8, height = 6)
#====================================================================================
# End of Script 
#====================================================================================


