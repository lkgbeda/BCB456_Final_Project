#!/bin/bash
#SBATCH --time=20:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=128G   # Maximum memory per node
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --account=geetu-lab
#SBATCH --job-name="sort"
#SBATCH --mail-user=lkgbeda@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL 
#SBATCH -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Coordinate_Sorted/sort_%J.out
#SBATCH -e /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Coordinate_Sorted/sort_%J.error

# Load required module
module purge
module load samtools

# Define directories
SAM_DIR="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Histat"
OUTPUT_DIR="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Coordinate_Sorted/"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Navigate to the SAM directory
cd ${SAM_DIR}

# Convert and sort SAM files to BAM
for f in *.sam
do
  e="${f%.sam}.bam"  # Change .sam to .bam
  echo "Processing ${f} → ${OUTPUT_DIR}/${e}"
  
  samtools sort -o ${OUTPUT_DIR}/${e} -O bam -T ${OUTPUT_DIR}/temp_${f%.*} -@ 16 ${f}
done

# Navigate to the output directory
cd ${OUTPUT_DIR}

# Index BAM files
for f in *.bam
do
  echo "Indexing ${f}"
  samtools index ${f}
done

echo "SAM to BAM conversion and indexing completed."


#Summary of what  the code does 
#The script is used to convert .sam files into sorted .bam files using samtools, and it also indexes those .bam files.
#SLURM directives manage job execution on a cluster.
#The samtools module is loaded to perform SAM/BAM operations.
#The script loops over SAM files, sorts them, and indexes the BAM files while storing them in a specified output directory.

