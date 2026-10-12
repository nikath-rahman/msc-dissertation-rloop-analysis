#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 00:30:00
#SBATCH --job-name=genomewide_rpe1_rianseq

# Genome-wide summary for the six RPE1 libraries: peak count and total peak
# coverage (bp). Etoposide does not cut at a defined recognition sequence, so
# no DSB site list is available for these samples.

DISS=${DISS:-$PWD}
PEAKS=${DISS}/peaks/rianseq
RESULTS=${DISS}/results/rianseq
SUMMARY=${RESULTS}/rpe1_genomewide_summary.txt

mkdir -p ${RESULTS}
printf 'Sample\tTotal_Peaks\tTotal_Peak_bp\n' > ${SUMMARY}

SAMPLES=(
  rpe1_WT_DMSO rpe1_WT_etoposide rpe1_WT_etoposide_RNaseH
  rpe1_upf1ko_DMSO rpe1_upf1ko_etoposide rpe1_upf1ko_etoposide_RNaseH
)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Counting peaks for ${SAMPLE}..."
    TOTAL=$(wc -l < ${PEAKS}/${SAMPLE}_peaks.narrowPeak)
    BP=$(awk '{sum+=$3-$2} END{print sum}' ${PEAKS}/${SAMPLE}_peaks.narrowPeak)
    printf '%s\t%s\t%s\n' "${SAMPLE}" "${TOTAL}" "${BP}" >> ${SUMMARY}
done

echo 'Genome-wide RPE1 peak comparison complete'
echo "Full summary written to ${SUMMARY}"
