#!/bin/bash
#SBATCH --time=20:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=128G   # Maximum memory per node
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --account=geetu-lab
#SBATCH --job-name="Index Sort"
#SBATCH --mail-user=lkgbeda@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL 
#SBATCH -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Index_Sorted/Indexsort_%J.out
#SBATCH -e /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Index_Sorted/Indexsort_%J.error

# Load required module
module purge
module load samtools

# Define directories
BAM_DIR="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Coordinate_Sorted/"
OUTPUT_DIR="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Index_Sorted/"

# Ensure output directory exists
mkdir -p ${OUTPUT_DIR}

# Navigate to the SAM directory
cd ${BAM_DIR}

# loop  through the BAM files
for f in *.bam
do
  e="${f%.bam}.bam"  # Prefixing the bam files
  echo "Processing ${f} → ${OUTPUT_DIR}/${e}"

# Coordinate Sorting
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

echo "BAM indexing completed."


#Summary of what  the code does 
#The script is used to convert .sam files into sorted .bam files using samtools, and it also indexes those .bam files.
#SLURM directives manage job execution on a cluster.
#The samtools module is loaded to perform SAM/BAM operations.
#The script loops over SAM files, sorts them, and indexes the BAM files while storing them in a specified output directory.

