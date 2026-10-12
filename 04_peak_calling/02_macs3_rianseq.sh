#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 02:00:00
#SBATCH --job-name=macs3_rianseq

# MACS3 peak calling for the 13 RIAN-seq libraries, with no input control.
# --extsize is the per-sample mean insert size from 03_insertsize_rianseq.sh.

DISS=${DISS:-$PWD}
FILT=${DISS}/filtered/rianseq
PEAKS=${DISS}/peaks/rianseq

module load Miniforge3/24.11.3-0
source ${EBROOTMINIFORGE3}/etc/profile.d/conda.sh
conda activate myenv
mkdir -p ${PEAKS}

# sample:extsize
SAMPLES=(
  "diva_siCtrl_DMSO:157"
  "diva_siCtrl_tam:112"
  "diva_siCtrl_DMSO_RNaseH:154"
  "diva_siCtrl_tam_RNaseH:117"
  "diva_siUpf1_DMSO:118"
  "diva_siUpf1_tam:145"
  "diva_siUpf1_tam_RNaseH:112"
  "rpe1_WT_DMSO:129"
  "rpe1_WT_etoposide:113"
  "rpe1_WT_etoposide_RNaseH:124"
  "rpe1_upf1ko_DMSO:130"
  "rpe1_upf1ko_etoposide:152"
  "rpe1_upf1ko_etoposide_RNaseH:139"
)

for ENTRY in "${SAMPLES[@]}"; do
    SAMPLE="${ENTRY%%:*}"
    EXTSIZE="${ENTRY##*:}"
    echo "Calling peaks for ${SAMPLE} (extsize=${EXTSIZE})..."
    macs3 callpeak \
        -t ${FILT}/${SAMPLE}.dedup.bam \
        -f BAMPE \
        -g hs \
        -n ${SAMPLE} \
        --outdir ${PEAKS} \
        --nomodel --extsize ${EXTSIZE} \
        -q 0.05
    echo "${SAMPLE} peak calling complete"
done

echo 'MACS3 peak calling complete for all 13 RIAN-seq samples (per-sample extsize)'
