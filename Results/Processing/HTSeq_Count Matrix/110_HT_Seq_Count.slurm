#!/bin/bash
#SBATCH --time=20:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=128G   # Maximum memory per node
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --account=geetu-lab
#SBATCH --job-name="HT-Seq"
#SBATCH --mail-user=lkgbeda@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL 
#SBATCH -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/HTSeq_Count/HTSeq_%J.out
#SBATCH -e /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/HTSeq_Count/HTSeq_%J.error

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module purge
module load python
source /work/LAS/geetu-lab/lkgbeda/PyEnv/bin/activate
module load py-htseq

# Define input directory and GTF file
gtf_file="//work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Ref_Genomes/Agrobacterium_fabrum/protein_coding_genes.gtf"
bam_dir="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Name_Sorted/"
out_dir="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/HTSeq_Count/"
mkdir -p "$out_dir"

# Loop through all BAM files and run htseq-count
for bam in "$bam_dir"/*.bam; do
  sample=$(basename "$bam" .bam)
  echo "Processing $sample..."
  htseq-count \
    --format=bam \
    --order=pos \
    --stranded=no \
    --type=gene \
    --idattr=locus_tag \
    --additional-attr=old_locus_tag \
    --mode=union \
    "$bam" "$gtf_file" > "$out_dir/${sample}_counts.txt"
done
