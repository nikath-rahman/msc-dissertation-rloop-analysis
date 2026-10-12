#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 02:00:00
#SBATCH --job-name=macs3

# MACS3 peak calling for spKAS-seq: enriched library against its input control.

DISS=${DISS:-$PWD}
ALIGNED=${DISS}/aligned
PEAKS=${DISS}/peaks

module load Miniforge3/24.11.3-0
source ${EBROOTMINIFORGE3}/etc/profile.d/conda.sh
conda activate myenv
mkdir -p ${PEAKS}

echo 'Running MACS3 peak calling...'
macs3 callpeak \
    -t ${ALIGNED}/Enrich_S2.filtered.bam \
    -c ${ALIGNED}/Input_S1.filtered.bam \
    -f BAMPE \
    -g 2948000000 \
    -n Enrich_S2_vs_Input_S1 \
    --outdir ${PEAKS} \
    -q 0.05 \
    --keep-dup auto

echo 'MACS3 peak calling complete'
