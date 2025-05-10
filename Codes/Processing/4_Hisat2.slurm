#!/bin/bash
#SBATCH --job-name="hisat"
#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=128G        # Maximum memory per node
#SBATCH --nodes=1         # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --account=geetu-lab
#SBATCH --mail-user=lkgbeda@iastate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Histat/hisat_%J.Out
#SBATCH -e /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Histat/hisat_%J.Error

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load hisat2
# Version: 2.1.0

# Set input and output directories
INPUT_DIR="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/TrimedSeq"
OUTPUT_DIR="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Histat"
GENOME_INDEX="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Ref_Genomes/Agrobacterium_fabrum/agrobacterium_fabrum_index"

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop over all *_1_paired.fastq.gz files and align as paired-end
for paired_fastq in "$INPUT_DIR"/*_1_paired.fastq.gz; do
    mate_fastq="${paired_fastq/_1_paired/_2_paired}"
    sample=$(basename "$paired_fastq" _1_paired.fastq.gz)

    echo "Running HISAT2 for sample: $sample"

    hisat2 -p 16 --rna-strandness R -x "$GENOME_INDEX" \
           -1 "$paired_fastq" -2 "$mate_fastq" \
           -S "$OUTPUT_DIR/${sample}.sam"

    echo "Alignment completed for: $sample"
done

echo "All paired-end alignments completed!"
