import csv
import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
print(str(Path(__file__).resolve().parent.parent.parent.parent))
from utils import ROOT_DIR
from utils.config import get_config

# Load chr_pos and ann
# Set input directory for JSON files and output file for BED format
config_obj = get_config()
wd = Path(config_obj["relevant_paths"]["generated_datasets"])
chr_pos = wd / '2024_10_22' / 'hg38positions.bed'  # Write a util to get latest dir
ann = '/data/shared/repos/biomuta-old/pipeline/convert_step2/liftover/Homo_sapiens.GRCh38.113.gff3'  # GFF file with genomic features
output_file = wd / '2024_10_22' / 'chr_pos_to_ensp.csv'

print("Starting process...")

# Read positions from chr_pos
positions = []
print("Loading SNP positions from h38positions.bed...")
with open(chr_pos, 'r') as f:
    for line in f:
        chrom, start, end = line.strip().split('\t')
        positions.append((chrom, int(start), int(end)))
print(f"Loaded {len(positions)} SNP positions.")

# Parse ann and group by chromosome
annotations = {}
print("Loading annotations from the gff file...")
with open(ann, 'r') as f:
    for line in f:
        if not line.startswith('#'):
            fields = line.strip().split('\t')
            chrom = 'chr' + fields[0]  # Add 'chr' prefix to match format in file1
            feature_start = int(fields[3])
            feature_end = int(fields[4])
            # Store annotations by chromosome
            if chrom not in annotations:
                annotations[chrom] = []
            annotations[chrom].append((feature_start, feature_end, line.strip()))
print(f"Loaded annotations for {len(annotations)} chromosomes.")

# Find matches and write to CSV
print("Finding matches and writing results to CSV...")
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(['chrom', 'start_pos', 'end_pos', 'ENSP'])  # Write header

    match_count = 0
    for chrom, snp_start, snp_end in positions:
        if chrom in annotations:
            for feature_start, feature_end, annotation in annotations[chrom]:
                if snp_start >= feature_start and snp_end <= feature_end:
                    # Extract ENSP from the annotation
                    ensp = None
                    for field in annotation.split(';'):
                        if 'protein_id=' in field:
                            ensp = field.split('=')[1]
                            break
                    
                    if ensp:
                        writer.writerow([chrom, snp_start, snp_end, ensp])
                        match_count += 1

    print(f"Finished writing to {output_file}. Total matches found: {match_count}")

print("Process completed.")