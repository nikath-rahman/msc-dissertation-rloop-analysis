#!/bin/bash

# Extracts the CHM13v2.0 coordinates (column 10) from the supplied DSB
# coordinate file and writes a BED file for BEDTools and regioneR.
# Ten of the 263 listed AsiSI sites have CHM13 coordinates.

DISS=${DISS:-$PWD}
INPUT=${DISS}/data/DSB_coordinates.txt
OUTPUT=${DISS}/data/DSB_coordinates_CHM13.bed

awk -F'\t' '
NR > 1 {
    coord = $10
    if (coord != "") {
        split(coord, a, ":")
        chrom = a[1]
        split(a[2], b, "-")
        print chrom "\t" b[1] "\t" b[2]
    }
}
' ${INPUT} > ${OUTPUT}

echo "Done. Output written to ${OUTPUT}"
wc -l ${OUTPUT}
head -5 ${OUTPUT}
