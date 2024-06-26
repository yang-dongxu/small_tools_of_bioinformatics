import os 
import argparse
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate an IGV batch script from a gene list and a GTF file.')
    parser.add_argument('--gene_list', required=True, help='Path to the gene list file.')
    parser.add_argument('--gtf_file', required=True, help='Path to the GTF file.')
    parser.add_argument('--flank', type=int, default=500, help='Flank width around the gene coordinates (in bp).')
    parser.add_argument('--output_script', default="./igv_batch.txt", help='Path to the output IGV batch script.')
    parser.add_argument('--session', required=True, help='Path to the IGV session file.')
    parser.add_argument('--snapshot_dir', default="./snapshot", help='Directory to save snapshots.')
    parser.add_argument('--exit', action="store_true", help='exit after batch command executed')

    return parser.parse_args()

def parse_gtf(gtf_file, genes):
    # Load GTF data using pandas
    col_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    gtf_data = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, names=col_names, dtype={'start': int, 'end': int})

    # Filter for transcripts
    transcripts = gtf_data[gtf_data['feature'] == 'transcript']

    # Extract gene_name and calculate transcript length
    transcripts['gene_name'] = transcripts['attributes'].str.extract('gene_name "([^"]+)"')
    transcripts['transcript_length'] = transcripts['end'] - transcripts['start']

    # Filter transcripts for genes of interest
    transcripts = transcripts[transcripts['gene_name'].isin(genes)]

    # Get the longest transcript for each gene
    idx = transcripts.groupby('gene_name')['transcript_length'].idxmax()
    longest_transcripts = transcripts.loc[idx]

    # Create dictionary of gene coordinates
    gene_coordinates = {row['gene_name']: (row['seqname'], row['start'], row['end'], row['strand'])
                        for index, row in longest_transcripts.iterrows()}

    return gene_coordinates
    
def write_igv_script(gene_coords, flank, output_script, session_path, snapshot_dir, exit):
    with open(output_script, 'w') as file:
        file.write('new\n')
        file.write(f'load {session_path}\n')
        file.write(f'snapshotDirectory {snapshot_dir}\n')
        for gene, coords in gene_coords.items():
            start = max(0, coords[1] - flank)
            end = coords[2] + flank
            file.write(f'goto {coords[0]}:{start}-{end}\n')
            file.write(f'snapshot {gene}.png\n')
        if exit:
            file.write("\nexit\n")
    
def main():
    args = parse_arguments()
    if not os.path.exists(args.snapshot_dir):
        os.makedirs(args.snapshot_dir)
    genes = set(line.strip() for line in open(args.gene_list))
    gene_coords = parse_gtf(args.gtf_file, genes)
    write_igv_script(gene_coords, args.flank, args.output_script, args.session, args.snapshot_dir, args.exit)
    print("IGV batch script has been generated successfully.")

if __name__ == "__main__":
    main()
