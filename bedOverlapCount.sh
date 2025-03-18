#!/bin/bash

# bedOverlapCount.sh - A tool to calculate overlap counts between multiple BED files
# Usage: ./bedOverlapCount.sh -r reference.bed -q "query_pattern.bed" -o output.tsv

# Default values
OUTFILE="overlap_counts.tsv"
HEADER=true
QUIET=false
KEEP_TEMP=false

# Help function
usage() {
    cat << EOF
bedOverlapCount - Calculate overlap counts between one reference and multiple query BED files

Usage: $(basename $0) [OPTIONS] -r REFERENCE.BED -q "QUERY_PATTERN" [QUERY2.BED ...]

Calculate overlap counts between multiple query BED files and a single reference BED file.
Outputs a TSV file with all columns from the reference file plus count columns for each query.

Options:
  -r, --reference FILE     Reference BED file (required)
  -q, --query PATTERNS     One or more query BED files or patterns (required)
                           Supports shell wildcards (e.g., "samples/*.bed")
  -o, --output FILE        Output file name (default: overlap_counts.tsv)
  -n, --no-header          Don't include a header in the output
  -k, --keep-temp          Keep temporary files for debugging
  -h, --help               Display this help and exit
  -v, --verbose            Print progress information

Example:
  $(basename $0) -r peaks.bed -q "sample*.bed" -o results.tsv
  $(basename $0) -r peaks.bed -q "data/exp1/*.bed" "data/exp2/*.bed" -o results.tsv
EOF
    exit 1
}

# Parse command line arguments
QUERY_PATTERNS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--reference)
            REFERENCE="$2"
            shift 2
            ;;
        -q|--query)
            shift
            while [[ $# -gt 0 && ! $1 =~ ^- ]]; do
                QUERY_PATTERNS+=("$1")
                shift
            done
            ;;
        -o|--output)
            OUTFILE="$2"
            shift 2
            ;;
        -n|--no-header)
            HEADER=false
            shift
            ;;
        -k|--keep-temp)
            KEEP_TEMP=true
            shift
            ;;
        -h|--help)
            usage
            ;;
        -v|--verbose)
            QUIET=false
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check for required arguments
if [ -z "$REFERENCE" ]; then
    echo "Error: Reference file is required."
    usage
fi

if [ ${#QUERY_PATTERNS[@]} -eq 0 ]; then
    echo "Error: At least one query file pattern is required."
    usage
fi

# Function to log messages if not quiet
log() {
    if [ "$QUIET" = false ]; then
        echo "$@"
    fi
}

# Check if reference file exists
if [ ! -f "$REFERENCE" ]; then
    echo "Error: Reference file '$REFERENCE' not found."
    exit 1
fi

# Expand query patterns into actual files
QUERIES=()
for PATTERN in "${QUERY_PATTERNS[@]}"; do
    # If pattern contains wildcard, expand it
    if [[ "$PATTERN" == *\** ]]; then
        # Use compgen to expand the pattern without word splitting
        MATCHES=$(compgen -G "$PATTERN" || echo "")
        if [ -z "$MATCHES" ]; then
            echo "Warning: Pattern '$PATTERN' didn't match any files."
            continue
        fi
        
        # Add each match to the queries array
        while IFS= read -r MATCH; do
            QUERIES+=("$MATCH")
        done <<< "$MATCHES"
    else
        # Not a pattern, just a regular file
        if [ ! -f "$PATTERN" ]; then
            echo "Error: Query file '$PATTERN' not found."
            exit 1
        fi
        QUERIES+=("$PATTERN")
    fi
done

# Check if we have any queries after expansion
if [ ${#QUERIES[@]} -eq 0 ]; then
    echo "Error: No valid query files found after expanding patterns."
    exit 1
fi

log "Found ${#QUERIES[@]} query files to process."

# Create a temporary directory
TEMP_DIR=$(mktemp -d)
trap 'if [ "$KEEP_TEMP" = false ]; then rm -rf "$TEMP_DIR"; fi' EXIT

# ---- FIX: Complete rewrite of the processing section ----

# First determine if reference has a header line
HAS_HEADER=false
if head -n 1 "$REFERENCE" | grep -q "^#"; then
    HAS_HEADER=true
fi

# Process each query file and get the counts
for i in "${!QUERIES[@]}"; do
    QUERY="${QUERIES[$i]}"
    BASENAME=$(basename "$QUERY" | sed 's/\.[^.]*$//')
    log "Processing query file: $QUERY"
    
    # Run bedtools intersect with counts
    bedtools intersect -a "$REFERENCE" -b "$QUERY" -c > "$TEMP_DIR/counts_${i}.bed"
    
    # If this is the first query, keep all columns
    if [ $i -eq 0 ]; then
        cp "$TEMP_DIR/counts_${i}.bed" "$TEMP_DIR/result.bed"
    else
        # Join with previous results using the reference columns as key
        # Count the number of columns in the reference file (minus the header if it exists)
        if [ "$HAS_HEADER" = true ]; then
            REF_COLS=$(head -n 2 "$REFERENCE" | tail -n 1 | awk '{print NF}')
        else
            REF_COLS=$(head -n 1 "$REFERENCE" | awk '{print NF}')
        fi
        
        # Extract just the count column from the current result
        cut -f$((REF_COLS+1)) "$TEMP_DIR/counts_${i}.bed" > "$TEMP_DIR/count_col_${i}.txt"
        
        # Paste this column to the existing result
        paste "$TEMP_DIR/result.bed" "$TEMP_DIR/count_col_${i}.txt" > "$TEMP_DIR/new_result.bed"
        mv "$TEMP_DIR/new_result.bed" "$TEMP_DIR/result.bed"
    fi
done

# Create the header for the output file
if [ "$HEADER" = true ]; then
    # Get column names from reference
    if [ "$HAS_HEADER" = true ]; then
        # Use the reference header line (without the #)
        HEAD_LINE=$(head -n 1 "$REFERENCE" | sed 's/^#//')
    else
        # Create generic column names for reference
        REF_COLS=$(head -n 1 "$REFERENCE" | awk '{print NF}')
        HEAD_LINE=""
        for ((i=1; i<=REF_COLS; i++)); do
            if [ -n "$HEAD_LINE" ]; then
                HEAD_LINE="${HEAD_LINE}\tCol${i}"
            else
                HEAD_LINE="Col${i}"
            fi
        done
    fi
    
    # Add query file names as column headers
    for QUERY in "${QUERIES[@]}"; do
        BASENAME=$(basename "$QUERY" | sed 's/\.[^.]*$//')
        HEAD_LINE="${HEAD_LINE}\t${BASENAME}"
    done
    
    # Write the header to the output file
    echo -e "$HEAD_LINE" > "$OUTFILE"
    
    # Write the data, skipping the header line if reference had one
    if [ "$HAS_HEADER" = true ]; then
        tail -n +2 "$TEMP_DIR/result.bed" >> "$OUTFILE"
    else
        cat "$TEMP_DIR/result.bed" >> "$OUTFILE"
    fi
else
    # No header, just copy the result
    cat "$TEMP_DIR/result.bed" > "$OUTFILE"
fi

log "Results written to $OUTFILE"
log "Processed ${#QUERIES[@]} query files against reference $REFERENCE"
