#!/bin/bash
#SBATCH --time=20:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=128G        # Maximum memory per node
#SBATCH --nodes=1         # number of nodes
#SBATCH --ntasks-per-node=16  # 16 processor core(s) per node
#SBATCH --account=geetu-lab
#SBATCH --job-name="Name_Sort"
#SBATCH --mail-user=lkgbeda@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL 
#SBATCH -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Name_Sorted/NameSorted_%J.out
#SBATCH -e /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Name_Sorted/NameSorted_%J.error

# LOAD MODULES
module purge
module load samtools

# Set working directory to where the BAM files are
cd /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Index_Sorted/

# Create output directory if it doesn't exist
mkdir -p /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Name_Sorted/

# Loop through all BAM files (excluding .bai index files)
for f in *.bam; do
    # Skip index files
    if [[ "$f" == *.bam.bai ]]; then
        continue
    fi

    # Set output filename
    out="${f%.bam}.bam"

    # Perform name sorting
    samtools sort -n -@ 16 -O bam -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Name_Sorted/$out \
        -T /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Name_Sorted/temp_$out \
        "$f"
done

echo "Name sorting complete for all BAM files!"
