#!/bin/bash

#SBATCH -J Filter_BAM
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5000M
#SBATCH --output=/your/path/ErrorOutput/FilterBAMlog.out
#SBATCH --error=/your/path/ErrorOutput/FilterBAMlog.err
#SBATCH --partition=....


# Define input CSV file
INPUT_FILE="/your/path/Data/samplesheetBAM.csv"
OUTPUT_BASE_DIR="/your/path/Data/Samples"

# Check if the input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: File '$INPUT_FILE' not found!"
    exit 1
fi

# Read the header and find the sample column (assuming 'SampleID' is the first column and 'dataDir' is the third column)
HEADER=$(head -n 1 "$INPUT_FILE")
echo "$HEADER"

while IFS=, read -r SAMPLE_ID DATA_DIR REST || [[ -n "$SAMPLE_ID" ]]; do
    # Skip the header row
    if [[ "$SAMPLE_ID" == "SampleID" ]]; then
        continue
    fi
    
    # Define output directory for this sample
    SAMPLE_OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE_ID}/MT_depth"
    mkdir -p "$SAMPLE_OUTPUT_DIR"
    
    # Construct the full path to the BAM file
    BAM_FILE="${DATA_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.bam"
    OUTPUT_FILE="${SAMPLE_OUTPUT_DIR}/${SAMPLE_ID}.sorted.bam"
    


    # Run samtools command
    if [[ -f "$BAM_FILE" ]]; then
        samtools view -h "$BAM_FILE" chrM > "$OUTPUT_FILE"

    else
        echo "Warning: BAM file for ${SAMPLE_ID}.bam not found in ${DATA_DIR}/${SAMPLE_ID}. Skipping."
    fi

done < "$INPUT_FILE"
