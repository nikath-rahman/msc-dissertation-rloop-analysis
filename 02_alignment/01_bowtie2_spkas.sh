#!/bin/bash
#SBATCH -p compute
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH -t 08:00:00
#SBATCH --job-name=bowtie2

# Builds the Bowtie2 index for CHM13v2.0 (if absent) and aligns the
# spKAS-seq libraries. Outputs coordinate-sorted BAM plus flagstat.

DISS=${DISS:-$PWD}
REF=${DISS}/reference/chm13v2.0
TRIM=${DISS}/trimmed
ALIGNED=${DISS}/aligned
LOGS=${DISS}/logs

module load Bowtie2/2.5.4-GCC-13.3.0
module load SAMtools/1.21-GCC-13.3.0
mkdir -p ${ALIGNED} ${LOGS}

if [ ! -f ${REF}.1.bt2 ] && [ ! -f ${REF}.1.bt2l ]; then
    echo 'Building Bowtie2 index from CHM13-T2T...'
    bowtie2-build --threads 8 ${REF}.fa ${REF}
else
    echo 'Index already exists. Skipping build step.'
fi

for SAMPLE in Enrich_S2 Input_S1; do
    echo "Aligning ${SAMPLE}..."
    bowtie2 -x ${REF} \
        -1 ${TRIM}/${SAMPLE}_merge_R1_val_1.fq.gz \
        -2 ${TRIM}/${SAMPLE}_merge_R2_val_2.fq.gz \
        -p 6 --no-unal --no-mixed --no-discordant -X 1000 \
        2> ${LOGS}/${SAMPLE}.bowtie2.log \
        | samtools sort -@ 2 -o ${ALIGNED}/${SAMPLE}.sorted.bam
    samtools index ${ALIGNED}/${SAMPLE}.sorted.bam
    samtools flagstat ${ALIGNED}/${SAMPLE}.sorted.bam > ${LOGS}/${SAMPLE}.flagstat.txt
    echo "${SAMPLE} alignment complete"
done

echo 'Bowtie2 pipeline complete'
