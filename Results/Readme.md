# Description of Folder Contents

## `Galaxy`

This folder contains files genrated from analysis using galaxy. Subfolder include `DEGs`, `Ensembl genome reference/ncbi_dataset`, `Figure 3a`, `WORKFLOW`, `csv files`, and `figure 4b`
`DEGs` : contains plots and diagrams illustrating differentially expressed genes
`Ensembl genome reference/ncbi_dataset` : contains the reference genome .gff file
`Figure 3a` : contains results of Expression Pattern Analysis of Genes Related to Environmental Information Processing
`WORKFLOW`: contains the pdf decribing the workflow of analysis
`csv files`: contains compiled dataset of upregulated and downregulated genes. It also contains the old and new nomencleture for expressed genes
`figure 4b` : contains results of the Expression Pattern Analysis of Genes Related to Cellular Processes

## `R`
This folder contains file output from r analysis
Subfolder include `DEG_Results` and `Figures`
`DEG_Results` : contains tab separated data of differentially expressed genes
`Figure` : contains all figures generated from R analysis. These include `down and upregulated genes diagrams`, `rldpca plots`,`FPKM heatmaps`, `FPKM distribution violin plots`, `FPKM distribution plots` and `gene exxpression summary`

## `Processing`
This folder contains intermediate files generated from data processing steps
Subfolders include `HTSeq_Count Matrix`, `Trimmed_Fastqc_Multiqc`, and `Untrimmed_Fastqc_Multiqc`
`HTSeq_Count Matrix` :contain count matrixes of gene expression
`Trimmed_Fastqc_Multiqc` : contains QC reports of read quality across genes

## `Paper Data`
This folder contains FPKM values from the paper. These values were analysed to confirm the reads or gene expressions reported in the paper.
