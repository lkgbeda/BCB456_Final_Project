#!/bin/bash
#SBATCH -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/TrimedSeq.Out 
#SBATCH -e /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/TrimedSeq.Errors
#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=128G        # Maximum memory per node
#SBATCH --nodes=1         # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --account=geetu-lab
#SBATCH --mail-user=lkgbeda@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load parallel
module load trimmomatic

# Set directories
WORK_DIR="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/"
OUTPUT_DIR="/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/TrimedSeq/"

# Create the output directory if it doesn’t exist
mkdir -p "$OUTPUT_DIR"

# Run Trimmomatic
cd "$WORK_DIR"

# Sample list from SRR15969054 to SRR15969071 (S1 to S18)
samples=(
SRR15969054
SRR15969055
SRR15969056
SRR15969057
SRR15969058
SRR15969059
SRR15969060
SRR15969061
SRR15969062
SRR15969063
SRR15969064
SRR15969065
SRR15969066
SRR15969067
SRR15969068
SRR15969069
SRR15969070
SRR15969071
)

# Loop through each sample and run Trimmomatic
for i in "${!samples[@]}"; do
    sample=${samples[$i]}
    sample_id="S$((i+1))"
    echo "# $sample_id - Paired-End Trimming"
    trimmomatic PE -phred33 \
    "${sample}_1.fastq.gz" "${sample}_2.fastq.gz" \
    "${OUTPUT_DIR}${sample}_1_paired.fastq.gz" "${OUTPUT_DIR}${sample}_1_unpaired.fastq.gz" \
    "${OUTPUT_DIR}${sample}_2_paired.fastq.gz" "${OUTPUT_DIR}${sample}_2_unpaired.fastq.gz" \
    ILLUMINACLIP:/work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/TruSeq3-PE-2.fa:2:30:10 \
    LEADING:28 TRAILING:30 MINLEN:36
done






