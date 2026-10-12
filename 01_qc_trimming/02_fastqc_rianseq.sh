#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 01:00:00
#SBATCH --job-name=01fastqc_rianseq

# FastQC on the raw RIAN-seq reads (13 libraries, GN_1 to GN_13).

DISS=${DISS:-$PWD}
RAWDIR=${DISS}/data/rianseq/X204SC26081427-Z01-F001/01.RawData
QC=${DISS}/qc/rianseq

module load FastQC/0.12.1-Java-11
mkdir -p ${QC}

FILES=()
for i in $(seq 1 13); do
    FILES+=( $(find ${RAWDIR}/GN_${i} -name '*_1.fq.gz') )
    FILES+=( $(find ${RAWDIR}/GN_${i} -name '*_2.fq.gz') )
done

fastqc "${FILES[@]}" --threads 4 -o ${QC}

echo 'FastQC complete - reports written to qc/rianseq/'
