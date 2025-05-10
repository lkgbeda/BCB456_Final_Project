#!/bin/bash
#SBATCH -o /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Fastqc.out
#SBATCH -e /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/Fastqc.errors
#SBATCH --time=10:00:00   # walltime limit (HH:MM:SS)
#SBATCH --mem=128G   # Maximum memory per node
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node
#SBATCH --account=geetu-lab
#SBATCH --mail-user=lkgbeda@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load fastqc
source /work/LAS/geetu-lab/lkgbeda/PyEnv/bin/activate
module load multiqc  # Ensure MultiQC is available


# Lane 1
cd /work/LAS/geetu-lab/lkgbeda/RNA_Seq/BCB546_Project/
fastqc *fastq.gz
multiqc . -o ../multiqc

