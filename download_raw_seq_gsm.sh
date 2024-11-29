#!/bin/bash

# Check if the necessary tools are installed
if ! command -v esearch &> /dev/null || ! command -v efetch &> /dev/null || ! command -v fasterq-dump &> /dev/null; then
    echo "Error: EDirect and SRA Toolkit must be installed and in your PATH."
    echo "for EDirect:  https://www.ncbi.nlm.nih.gov/books/NBK179288/"
    echo "for SRA Toolkit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit" 
    exit 1
fi

# Default values for parameters
threads=16
output_dir="RawData"

# Function to print usage information
usage() {
    echo "Usage: $0 -i <gsm_file> [-o <output_dir>] [-t <threads>]"
    echo "  -i <gsm_file>     : File containing list of GSM IDs, one per line."
    echo "  -o <output_dir>   : Directory to save FASTQ files (default: fastq_files)."
    echo "  -t <threads>      : Number of threads for fasterq-dump (default: 4)."
    echo "  -h                : Display this help message."
    exit 1
}


# Parse command-line arguments
while getopts "i:o:t:h" opt; do
    case "${opt}" in
        i)
            gsm_file="${OPTARG}"
            ;;
        o)
            output_dir="${OPTARG}"
            ;;
        t)
            threads="${OPTARG}"
            ;;
        h)
            usage
            ;;
        *)
            usage
            ;;
    esac
done
gsm_file=$(realpath "$gsm_file")

# Check if the GSM file is specified
if [ -z "${gsm_file}" ]; then
    echo "Error: GSM file is not specified."
    usage
fi

# Create and navigate to the output directory
mkdir -p "$output_dir"
cd "$output_dir" || exit

# Loop through each GSM ID and process
while read -r GSM_ID <&3; do
    echo "Processing $GSM_ID..."
    
    # Fetch corresponding SRX accession numbers
    SRX_IDS=$(esearch -db gds -query "$GSM_ID" | \
              efetch -format docsum | \
              xtract -pattern ExtRelation -block ExtRelation -element TargetObject | \
              grep '^SRX')
    
    if [[ -z "$SRX_IDS" ]]; then
        echo "Warning: No SRX IDs found for $GSM_ID"
        continue
    fi

    if [[ $(echo "$SRX_IDS" | wc -l) -gt 1 ]]; then
        echo "Warning: Multiple SRX IDs ($SRX_IDS) found for $GSM_ID."
    fi
    
    # For each SRX ID, find the corresponding SRR IDs
    for SRX_ID in $SRX_IDS; do
        echo "Fetching SRR IDs for $SRX_ID..."
        
        SRR_IDS=$(esearch -db sra -query "$SRX_ID" | \
                  efetch -format runinfo | \
                  cut -d ',' -f 1 | grep '^SRR')
                
        echo "SRR IDs for $GSM_ID : $SRR_IDS"
        
        if [[ -z "$SRR_IDS" ]]; then
            echo "Warning: No SRR IDs found for $SRX_ID"
            continue
        fi
        
        # Download FASTQ files for each SRR ID using fasterq-dump
        for SRR_ID in $SRR_IDS; do
            echo "Downloading FASTQ for $GSM_ID $SRX_ID $SRR_ID using fasterq-dump..."
            fasterq-dump --split-files -m 1024MB -p  -e "$threads" "$SRR_ID" -O $GSM_ID
        done
    done

done 3< "$gsm_file"

echo "All downloads complete."
