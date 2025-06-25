#!/bin/bash
#SBATCH -J MTGATK
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=2000M
#SBATCH --output=/your/path/ErrorOutput/gatklog.out
#SBATCH --error=/your/path/ErrorOutput/gatklog.err
#SBATCH --partition=....

EmbryoID=$1

cd /your/path/Samples/$EmbryoID/MT_depth

module load bioinf/samtools/1.15.1


samtools index $EmbryoID.bam

module load bioinf/gatk/4.3.0.0

# Run GATK HaplotypeCaller
gatk --java-options "-Xmx8g" HaplotypeCaller \
    -R /your/path/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.pgt.fna \
    -I /your/path/Data/Samples/$EmbryoID/MT_depth/$EmbryoID.bam \
    -O /your/path/Data/Samples/$EmbryoID/MT_depth/$EmbryoID.g.vcf \
    -ERC GVCF \
    -L chrM
