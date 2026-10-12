#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 02:00:00
#SBATCH --job-name=trimgalore

# Adapter and quality trimming of the spKAS-seq libraries.
# Default Illumina adapter detection; Q20; minimum post-trim length 20 bp.

DISS=${DISS:-$PWD}
RAW=${DISS}/data/raw
TRIM=${DISS}/trimmed

module load Trim_Galore/0.6.10-GCCcore-12.3.0
module load cutadapt/4.9-GCCcore-12.3.0
mkdir -p ${TRIM}

for SAMPLE in Enrich_S2 Input_S1; do
    echo "Trimming ${SAMPLE}..."
    trim_galore \
        --paired \
        --quality 20 \
        --length 20 \
        --cores 4 \
        --fastqc \
        --output_dir ${TRIM} \
        ${RAW}/${SAMPLE}_merge_R1.fastq.gz \
        ${RAW}/${SAMPLE}_merge_R2.fastq.gz
done

echo 'Trim Galore complete'
