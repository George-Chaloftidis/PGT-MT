#!/bin/bash
#SBATCH -J MTdepth_sampleID
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=5000M
#SBATCH --output=/your/path/ErrorOutput/MTdepthlog.out
#SBATCH --error=/your/path/ErrorOutput/MTdepthlog.err
#SBATCH --partition=....

# Load modules 
module load bioinf/samtools

# Define the base directory and the input CSV file name
BASE_DIR="/your/path/Data"
INPUT_FILE="${BASE_DIR}/samplesheet.csv"
OUTPUT_BASE_DIR="${BASE_DIR}/Samples"

# Check if the input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: File '$INPUT_FILE' not found!"
    exit 1
fi

# Read the header and find the sample column (assuming 'SampleID' is the first column and 'dataDir' is the second column)
HEADER=$(head -n 1 "$INPUT_FILE")
echo "$HEADER"

# Loop through the CSV file
while IFS=, read -r SAMPLE_ID DATA_DIR REST || [[ -n "$SAMPLE_ID" ]]; do
    # Skip the header row
    if [[ "$SAMPLE_ID" == "SampleID" ]]; then
        continue
    fi

    # Define output directory for this sample
    SAMPLE_OUTPUT_DIR="${OUTPUT_BASE_DIR}/${SAMPLE_ID}/MT_depth"
    mkdir -p "$SAMPLE_OUTPUT_DIR"
    
    # Construct the full path to the BAM file
    BAM_FILE="${DATA_DIR}/${SAMPLE_ID}/MT_depth/${SAMPLE_ID}.bam"

    # Print the BAM file path to debug
    echo "Looking for BAM file at: ${BAM_FILE}"
    
    # Define the output file without the full path
    OUTPUT_FILE="${SAMPLE_ID}_MT.txt"

    # Check if the BAM file exists
    if [[ -f "$BAM_FILE" ]]; then
        echo "Found BAM file: $BAM_FILE"
        
        # Run samtools depth and modify the output format to only include the file name
        samtools depth -a -b chrM.bed -H "${BAM_FILE}" | \
        sed "s|/your/path/Data/Samples/${SAMPLE_ID}/MT_depth/${SAMPLE_ID}.bam|${SAMPLE_ID}.bam|" > "${SAMPLE_OUTPUT_DIR}/${SAMPLE_ID}_MT.txt"
    else
        echo "Warning: BAM file for ${SAMPLE_ID} not found in ${DATA_DIR}. Skipping."
    fi

done < "$INPUT_FILE"
