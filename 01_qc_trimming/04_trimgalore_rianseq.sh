#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 06:00:00
#SBATCH --job-name=trimgalore_rianseq

# Adapter and quality trimming of the 13 RIAN-seq libraries.
# --nextera because the libraries are Tn5-tagmented; Q20; minimum length 20 bp.
# GN number to sample name: see sample_sheet.tsv in the repository root.

DISS=${DISS:-$PWD}
RAWDIR=${DISS}/data/rianseq/X204SC26081427-Z01-F001/01.RawData
TRIM=${DISS}/trimmed/rianseq

module load Trim_Galore/0.6.10-GCCcore-12.3.0
module load cutadapt/4.9-GCCcore-12.3.0
mkdir -p ${TRIM}

for i in $(seq 1 13); do
    echo "Trimming GN_${i}..."
    trim_galore \
        --paired \
        --nextera \
        --quality 20 \
        --length 20 \
        --cores 4 \
        --fastqc \
        --output_dir ${TRIM} \
        $(find ${RAWDIR}/GN_${i} -name '*_1.fq.gz') \
        $(find ${RAWDIR}/GN_${i} -name '*_2.fq.gz')
done

echo 'Trim Galore complete for all 13 RIAN-seq samples'
