#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 01:00:00
#SBATCH --job-name=01fastqc

# FastQC on the raw spKAS-seq reads (enriched and input libraries).
# Set DISS to the project root before running.

DISS=${DISS:-$PWD}
RAW=${DISS}/data/raw
QC=${DISS}/qc

module load FastQC/0.12.1-Java-11
mkdir -p ${QC}

fastqc \
    ${RAW}/Enrich_S2_merge_R1.fastq.gz \
    ${RAW}/Enrich_S2_merge_R2.fastq.gz \
    ${RAW}/Input_S1_merge_R1.fastq.gz \
    ${RAW}/Input_S1_merge_R2.fastq.gz \
    --threads 4 \
    -o ${QC}

echo 'FastQC complete - reports written to qc/'
