#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH -t 01:00:00
#SBATCH --job-name=bedtools_diva_rianseq

# Peaks overlapping the ten DSB sites, for the seven DIvA RIAN-seq libraries.
# The DSB coordinates are AsiSI sites in DIvA cells, so RPE1 is not tested here.

DISS=${DISS:-$PWD}
PEAKS=${DISS}/peaks/rianseq
DSB=${DISS}/data/DSB_coordinates_CHM13.bed
RESULTS=${DISS}/results/rianseq
SUMMARY=${RESULTS}/diva_dsb_overlap_summary.txt

module load BEDTools/2.31.1-GCC-13.3.0
mkdir -p ${RESULTS}
printf 'Sample\tTotal_Peaks\tPeaks_Overlapping_DSB\n' > ${SUMMARY}

SAMPLES=(
  diva_siCtrl_DMSO diva_siCtrl_tam diva_siCtrl_DMSO_RNaseH diva_siCtrl_tam_RNaseH
  diva_siUpf1_DMSO diva_siUpf1_tam diva_siUpf1_tam_RNaseH
)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Intersecting ${SAMPLE} with DSB sites..."
    bedtools intersect -u \
        -a ${PEAKS}/${SAMPLE}_peaks.narrowPeak \
        -b ${DSB} > ${RESULTS}/${SAMPLE}_dsb_overlap.bed
    TOTAL=$(wc -l < ${PEAKS}/${SAMPLE}_peaks.narrowPeak)
    OVERLAP=$(wc -l < ${RESULTS}/${SAMPLE}_dsb_overlap.bed)
    printf '%s\t%s\t%s\n' "${SAMPLE}" "${TOTAL}" "${OVERLAP}" >> ${SUMMARY}
done

echo "BEDTools DSB overlap analysis complete for all DIvA samples"
echo "Full summary written to ${SUMMARY}"
