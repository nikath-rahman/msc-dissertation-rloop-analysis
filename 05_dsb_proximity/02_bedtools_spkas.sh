#!/bin/bash

# Overlap of spKAS-seq peaks with the ten DSB sites:
#  - exact overlap (bedtools intersect)
#  - peaks within 1 kb (bedtools window -w 1000)

DISS=${DISS:-$PWD}
PEAKS=${DISS}/peaks/Enrich_S2_vs_Input_S1_peaks.narrowPeak
DSB=${DISS}/data/DSB_coordinates_CHM13.bed
RESULTS=${DISS}/results

module load BEDTools/2.31.1-GCC-13.3.0
mkdir -p ${RESULTS}

echo "=== Exact overlap: peaks directly intersecting a DSB site ==="
bedtools intersect -a ${PEAKS} -b ${DSB} -wa -wb > ${RESULTS}/exact_overlap.bed
wc -l ${RESULTS}/exact_overlap.bed

echo "=== Window overlap: peaks within 1000 bp of a DSB site ==="
bedtools window -a ${PEAKS} -b ${DSB} -w 1000 > ${RESULTS}/window_overlap.bed
wc -l ${RESULTS}/window_overlap.bed

echo "=== DSB sites with a peak within 1000 bp ==="
bedtools window -a ${DSB} -b ${PEAKS} -w 1000 | cut -f1-3 | sort -u | wc -l

echo "Done. Results in ${RESULTS}/"
