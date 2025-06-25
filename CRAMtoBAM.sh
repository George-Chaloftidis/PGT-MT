#!/bin/bash

#SBATCH -J CRAM_to_BAM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=5000M
#SBATCH --output=/your/path/ErrorOutput/CRAMtoBAMlog.out
#SBATCH --error=/your/path/ErrorOutput/CRAMtoBAMlog.err
#SBATCH --partition=.....
#SBATCH --exclude=.....
#SBATCH --array=0-50  # Array range, adjust by number of CRAM files


# Conda environment activation
source activate /your/conda/directory

# Set Directory and Reference
DATA_DIR="your/path/Data/Samples"
REF_DIR="your/path/GRCh38_full_analysis_set_plus_decoy_hla.pgt.fa"

# List of CRAM files
samples=($(find "$DATA_DIR" -type f -name "*.cram"))

# Select CRAM file based on SLURM array task ID
SAMPLE="${samples[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$SAMPLE" .cram)

# Set output path
CRAM_DIR=$(dirname "$SAMPLE")
BAM="${CRAM_DIR}/${BASENAME}.bam"

# Convert CRAM to BAM
echo "Converting $SAMPLE to BAM..."
samtools view -b -@ 16 -T "${REF_DIR}" -o "$BAM" "$SAMPLE"

# Index the BAM
echo "Indexing $BAM..."
samtools index "$BAM"
