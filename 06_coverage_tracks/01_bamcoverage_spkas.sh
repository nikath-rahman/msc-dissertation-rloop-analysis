#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH -t 02:30:00
#SBATCH --job-name=deeptools

# RPGC-normalised bigWig tracks for the spKAS-seq libraries (used for the
# coverage figure at the DSB sites).

DISS=${DISS:-$PWD}
ALIGNED=${DISS}/aligned
BIGWIGS=${DISS}/bigwigs
GENOME_SIZE=2948000000

module load Miniforge3/24.11.3-0
source ${EBROOTMINIFORGE3}/etc/profile.d/conda.sh
conda activate myenv
mkdir -p ${BIGWIGS}

for SAMPLE in Enrich_S2 Input_S1; do
    echo "Generating bigWig for ${SAMPLE}..."
    bamCoverage \
        -b ${ALIGNED}/${SAMPLE}.filtered.bam \
        -o ${BIGWIGS}/${SAMPLE}.bw \
        --normalizeUsing RPGC \
        --effectiveGenomeSize ${GENOME_SIZE} \
        -p 4
done

echo 'bigWig generation complete'
