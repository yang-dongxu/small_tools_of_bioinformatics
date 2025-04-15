#!/bin/bash

# Check if the necessary tools are installed
if ! command -v esearch &> /dev/null || ! command -v efetch &> /dev/null ; then
    echo "Error: EDirect and SRA Toolkit must be installed and in your PATH."
    echo "For EDirect:  https://www.ncbi.nlm.nih.gov/books/NBK179288/"
    echo "For SRA Toolkit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit" 
    exit 1
fi

# Default values for parameters
threads=16
output_dir="RawData"
download_method="wget"  # Default download method
aspera_key="~/.aspera/connect/etc/asperaweb_id_dsa.openssh"
# Function to print usage information
usage() {
    echo "Usage: $0 -i <gsm_file> [-o <output_dir>] [-t <threads>] [-m <method>]"
    echo "  -i <gsm_file>     : File containing list of GSM IDs, one per line."
    echo "  -o <output_dir>   : Directory to save FASTQ files (default: RawData)."
    echo "  -m <method>       : Download method ('sra' or 'wget', or 'ascp', default: wget)."
    echo "  -t <threads>      : Number of threads for fasterq-dump (default: 16)."
    echo "  -k <ascp_key>     : Path to Aspera key file. Default: ~/.aspera/connect/etc/asperaweb_id_dsa.openssh"
    echo "  -h                : Display this help message."
    exit 1
}

# Parse command-line arguments
while getopts "i:o:t:m:k:h" opt; do
    case "${opt}" in
        i)
            gsm_file="${OPTARG}"
            ;;
        o)
            output_dir="${OPTARG}"
            ;;
        t)
            threads="${OPTARG}"
            if [[ ! "$threads" =~ ^[1-9][0-9]*$ ]]; then
                echo "Error: Number of threads must be a positive integer."
                exit 1
            fi
            ;;
        m)
            download_method="${OPTARG}"
            if [[ "$download_method" != "sra" && "$download_method" != "wget" && "$download_method" != "ascp" ]]; then
                echo "Error: Invalid method. Choose 'sra', 'wget' or 'ascp'. Choosed: '$download_method'. "
                exit 1
            fi
            ;;
        k)
            aspera_key="${OPTARG}"
            ;;
        h)
            usage
            ;;
        *)
            usage
            ;;
    esac
done


if  ! command -v fasterq-dump &> /dev/null ; then
    echo "Error: SRA Toolkit must be installed and in your PATH."
    echo "For SRA Toolkit:  https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit" 
    exit 1
fi

# Normalize GSM file's path
gsm_file=$(realpath "$gsm_file")

# Check if the GSM file is specified
if [ -z "${gsm_file}" ]; then
    echo "Error: GSM file is not specified."
    usage
fi

if [[ "$download_method" == "ascp" ]]; then
    if ! command -v ascp &> /dev/null ; then
        echo "Error: Aspera CLI (ascp) must be installed and in your PATH."
        echo "For Aspera CLI:  https://downloads.asperasoft.com/"
        exit 1
    fi
    if [[ -z "$aspera_key" ]]; then
        echo "Error: Aspera key file is not specified."
        echo "Use -k <ascp_key> to specify the path to the Aspera key file."
        exit 1
    fi
fi

# Create and navigate to the output directory
mkdir -p "$output_dir"
cd "$output_dir" || exit

# Loop through each line in the GSM file and process the first column
while read -r line <&3; do
    GSM_ID=$(echo "$line" | awk '{print $1}')
    echo "Processing $GSM_ID..."
    
    # Fetch corresponding SRX accession numbers
    SRX_IDS=$(esearch -db gds -query "$GSM_ID" | efetch -format docsum | xtract -pattern ExtRelation -block ExtRelation -element TargetObject | grep '^SRX')
    
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
        
        SRR_IDS=$(esearch -db sra -query "$SRX_ID" | efetch -format runinfo | cut -d ',' -f 1 | grep '^SRR')
        
        if [[ -z "$SRR_IDS" ]]; then
            echo "Warning: No SRR IDs found for $SRX_ID"
            continue
        fi
        
        echo "SRR IDs for $GSM_ID: $SRR_IDS"
        
        # Create a directory for the current GSM ID
        mkdir -p "$GSM_ID"
        
        # Download FASTQ files for each SRR ID
        if [[ "$download_method" == "sra" ]]; then
            for SRR_ID in $SRR_IDS; do
                echo "Downloading FASTQ for $GSM_ID using fasterq-dump for SRR ID: $SRR_ID..."
                fasterq-dump --split-files -p -e "$threads" "$SRR_ID" -O "$GSM_ID"
            done
        else
            for SRR_ID in $SRR_IDS; do
                url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${SRR_ID}&result=read_run&fields=study_accession,sample_accession,experiment_accession,run_accession,tax_id,scientific_name,library_strategy,fastq_bytes,fastq_ftp,fastq_aspera,submitted_ftp,sra_ftp,sample_title&format=tsv&download=true&limit=0"
                tsv_file="${GSM_ID}/${SRR_ID}.tsv"
                
                echo "Fetching: $url"
                if ! wget -q -c "$url" -O "$tsv_file"; then
                    echo "Error: Failed to download $tsv_file"
                    continue
                fi
                
                if [[ "$download_method" == "ascp" ]]; then
                    fastqs=$(cut -f 10 "$tsv_file" | tail -n +2 | tr ';' $'\n' ) # Skip 
                    echo "Found FASTQ URLs for $SRR_ID:"
                    for fastq in $fastqs; do
                        echo "Downloading $fastq using ascp..."
                        ascp -i "$aspera_key" -k 1 -T -l 300m -P33001 "era-fasp@${fastq}" "$GSM_ID"
                        if [[ $? -ne 0 ]]; then
                            echo "Error: Failed to download $fastq"
                        fi
                    done
                elif [[ "$download_method" == "wget" ]]; then
                    fastqs=$(cut -f 9 "$tsv_file" | tail -n +2 | tr ';' $'\n' ) # Skip header
                    echo "Found FASTQ URLs for $SRR_ID:"
                    
                    for fastq in $fastqs; do
                        if ! wget  -c "$fastq" -P "$GSM_ID"; then
                            echo "Error: Failed to download $fastq"
                        fi
                    done
                else
                    echo "Error: Invalid download method. Choose 'sra', 'wget' or 'ascp'."
                    exit 1
                fi
            done
        fi
    done

done 3< "$gsm_file"

echo "All downloads complete."
