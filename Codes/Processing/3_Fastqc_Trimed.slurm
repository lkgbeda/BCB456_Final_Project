#!/bin/bash
#SBATCH -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/TrimedSeq/Fastqc.out
#SBATCH -e /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/TrimedSeq/Fastqc.errors
#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=128G        # Maximum memory per node
#SBATCH --nodes=1         # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --account=geetu-lab
#SBATCH --mail-user=lkgbeda@iastate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load fastqc
source /work/LAS/geetu-lab/lkgbeda/PyEnv/bin/activate
module load multiqc

# Set directories
INPUT_DIR="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/TrimedSeq"
FASTQC_DIR="$INPUT_DIR/Fastqc"
MULTIQC_DIR="$FASTQC_DIR/multiqc"

# Create output directories
mkdir -p "$FASTQC_DIR"
mkdir -p "$MULTIQC_DIR"

# Run FastQC only on *_paired.fastq.gz files
fastqc "$INPUT_DIR"/*_paired.fastq.gz -o "$FASTQC_DIR"

# Run MultiQC on the FastQC results
multiqc "$FASTQC_DIR" -o "$MULTIQC_DIR"
