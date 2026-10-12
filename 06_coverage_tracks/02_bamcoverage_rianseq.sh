#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 02:00:00
#SBATCH --job-name=bamcoverage_rianseq

# CPM-normalised bigWig tracks (10 bp bins) for the 13 RIAN-seq libraries.

DISS=${DISS:-$PWD}
FILT=${DISS}/filtered/rianseq
BIGWIG=${DISS}/bigwig/rianseq

module load Miniforge3/24.11.3-0
source ${EBROOTMINIFORGE3}/etc/profile.d/conda.sh
conda activate myenv
mkdir -p ${BIGWIG}

SAMPLES=(
  diva_siCtrl_DMSO diva_siCtrl_tam diva_siCtrl_DMSO_RNaseH diva_siCtrl_tam_RNaseH
  diva_siUpf1_DMSO diva_siUpf1_tam diva_siUpf1_tam_RNaseH
  rpe1_WT_DMSO rpe1_WT_etoposide rpe1_WT_etoposide_RNaseH
  rpe1_upf1ko_DMSO rpe1_upf1ko_etoposide rpe1_upf1ko_etoposide_RNaseH
)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Generating CPM bigwig for ${SAMPLE}..."
    bamCoverage \
        -b ${FILT}/${SAMPLE}.dedup.bam \
        -o ${BIGWIG}/${SAMPLE}.cpm.bw \
        --normalizeUsing CPM \
        --binSize 10 \
        -p 4
done

echo 'bamCoverage normalisation complete for all 13 RIAN-seq samples'
