### 2. Author (s)- YEAR.md 
# Replication Report: Comparative Transcriptome Analysis of Agrobacterium tumefaciens Reveals the Molecular Basis for the Recalcitrant Genetic Transformation of Camellia sinesis L
Original Authors:  Ke Jin 1,2,†, Na Tian 1,2,†, Jorge Freire da Silva Ferreira 3, Devinder Sandhu 3, Lizheng Xiao 1, Meiyi Gu 1, Yiping Luo 2, Xiangqin Zhang 2, Guizhi Liu 2, Zhonghua Liu 1,2,*, Jianan Huang 1,2,*, Shuoqian Liu 1,2,*
Original Year: 2022 May, 11. 
Replication Authors: Lakpa Sherpa, Faizo Kasule, Kofi Antwi Appiagyei, Lucky Kofi Gbeda, and Rishaniya Parthasarathy 
Replication date: 09May2025
## Original Study 
The study studied why tea plat (Camellia sinesis) are resistant to Agrobacterium-mediated transformation (AMT), a critical tool for genetic modification. By comparing AMT in tea leaves to the highly transformable tobacco (Nicotiana benthamiana), the authors identify physical and molecular barrier in tea that disrupt Agrobacterium infection, including 
Method: 
1.	Bacterial culture: Agrobacterium tumefacien strain GV3101 was cultured in LB medium at 28C overnight. 
2.	Tea and tobacco leaf disc (0.5cm x 0.5cm) were soaked in bacterial suspension for 20 mis 
3.	Leaved discs were fixed at intervals (D0-D4) for SEM imaging via Hitachi SU8100 SEM. 
4.	Transcriptomics: Bacterial RNA was extracted from tea/tobacco co-cultures at D0, D3, D4. rRNA depleted libraries sequenced on Illumina Novaseq. 
Alignment to A. tumefaciens genome+ plasmid (Bowtie2). DEGs inedited (DESeq2, P-adv<0.05) and annotated (GO/KEGG). 
5.	qRT-PCR validation: RNA extracted reversed transcribed to cDNA amplified using TB Green Premix. It was normalized to housekeeping genes (gyrB, dnaC, atu8171). Correlation performed between RNA-seq and qRT-PCTR data (Pearson). 

## Replication Details 
Data analysis was done using Galaxy and confirmed using R. The paper provided only the BioProject number on NCBI. No count matrix or GitHub repository were available. 
Data from the paper that was replicated include:
•	Gene Ontology enrichment analysis
•	Differential gene expression
•	Expression Pattern Analysis of Genes related to Environmental Information Processing
•	Expression Pattern Analysis of genes related to Cellular Processes

### Data 
-	Source: Samples and references were downloaded from NCBI ENA where paper provided the BioProject number PRJNA764576 and reference genome “A. tumefaciens str. C58”
-	Change made: Any modification or cleaning steps 
o	FastQC & MultiQC 
o	Trimmed low-quality reads and adapters with Trimmomatic
-	
### Methods 
-	Workflow of RNA-Seq were created using HPC and R software.  
o	After quality control checks and trimming
o	Reads were aligned to reference genome using HISAT2
o	SAM files were converted to BAM files using Samtools
o	These were sorted by genomic coordinates followed by indexing
o	We continued to sort by name
o	HTSeq0count was used to quantify gene expression or count matrix
o	R was used to calculate FPKM
o	Normalization & PCA was done using R
o	DSeq2 was used to analyse differentially expressed genes in R
o	DEGs were filtered and visualized in R
-	Difference from original: any methodological variation
o	
o	The major variation was the kind of tools used in analysis
o	We used Galaxy and R
o	The paper 
-	Challenges encountered: Description of any obstacles 
## Results 
Summary of how results compare to the original findings 
-	Tables/figures successfully reproduced 
-	Any difference found 
-	Possible reason of discrepancies 
## Conclusion 
-	Similar DEG pattern were observed on day 3 and day 4, but not on day 0. 
-	Key genes, especially those involved in chemotaxis, quorum sensing, and flagellar assembly, were consistently downregulated in Agrobacterium exposed to tea.
-	The original paper could be strengthened by analyzing regulatory changes between time points to better understand dynamic gene expression over the infection timeline
-	Time course analysis
-	Paper didn’t provide read me file and code for reproducibility
-	Proper documentation lacking

