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

# Navigate to merged files directory
cd /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Name_Sorted/

# Loop through files and perform htseq-count

for i in *.bam
do
e=`echo ${i} | sed 's/.bam/_count.txt/g'`
echo ${e}

htseq-count -f bam -s reverse -m intersection-nonempty --secondary-alignments ignore --supplementary-alignments ignore  ${i} /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Ref_Genomes/Agrobacterium_fabrum/protein_coding_genes.gtf > /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/HTSeq_Count/${e}
done
