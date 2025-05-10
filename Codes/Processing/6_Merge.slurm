#!/bin/bash
#SBATCH --job-name="merge replicates"
#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=128G   # Maximum memory per node
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --account=geetu-lab
#SBATCH --mail-user=lkgbeda@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Merged/Merge_%J.Out
#SBATCH -e /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Merged/Merge_%J.Error

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# Load required module
module purge
module load samtools

cd /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Unmerged/

# Replace below IDs with your own sample IDs or adjust for loop so that you are merging the technical replicates. 
manual_groups=("CD0" "CD3" "CD4" "ND0" "ND3" "ND4")

for group in "${manual_groups[@]}"; do
    files=$(ls *_${group}.bam 2>/dev/null | grep -v ".bam.bai")
    if [ -n "$files" ]; then
        echo "Merging manually listed group: $group"
        samtools merge /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Merged/merged_${group}.bam $files
    fi
done
