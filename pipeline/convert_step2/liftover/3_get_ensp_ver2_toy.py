import csv
import pickle
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
from utils import ROOT_DIR
from utils.config import get_config

# Load input_bed and ann
# Set input directory for JSON files and output file for BED format
config_obj = get_config()
wd = Path(config_obj["relevant_paths"]["generated_datasets"])
input_bed = wd / '2024_10_22' / 'liftover' / 'hg38_combined_toy.bed'  # Write a util to get latest dir
ann_dir = Path(config_obj["relevant_paths"]["downloads"])
ann = ann_dir / 'ensembl' / 'Homo_sapiens.GRCh38.113.gff3'  # GFF file with genomic features
output_file = wd / '2024_10_22' / 'mapping_ids' / 'chr_pos_to_ensp_toy.csv'


# Step 1: Load and Serialize 'ann' Annotations
def parse_and_serialize_ann():
    annotations = {}
    print("Parsing annotations from ann...")
    with open(ann, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                fields = line.strip().split('\t')
                chrom = 'chr' + fields[0]  # Add 'chr' prefix to match format in input_bed
                feature_start = int(fields[3])
                feature_end = int(fields[4])
                if chrom not in annotations:
                    annotations[chrom] = []
                annotations[chrom].append((feature_start, feature_end, line.strip()))
    with open('annotations.pkl', 'wb') as f:
        pickle.dump(annotations, f)
    print("Serialized annotations to 'annotations.pkl'")

# Step 2: Load 'ann' from Serialized File
def load_annotations():
    with open('annotations.pkl', 'rb') as f:
        return pickle.load(f)

# Step 3: Process 'input_bed' and Write Results in Batches
def process_large_input_bed():
    # Load serialized annotations
    annotations = load_annotations()

    # Open output CSV file
    with open(output_file, 'w', newline='') as csvfile:
        # Define the new headers for the output file
        writer = csv.writer(csvfile)
        writer.writerow(['chr_id', 'start_pos', 'end_pos', 'entrez_gene_id', 'prot_change', 'ENSP'])

        batch = []
        batch_size = 10000  # Define batch size for writing

        print("Processing SNP positions from input_bed and writing to CSV...")

        with open(input_bed, 'r') as f:
            # Skip the header line (if the first line is a header)
            header_skipped = False

            for i, line in enumerate(f, start=1):
                # Skip the header line
                if not header_skipped:
                    header_skipped = True
                    continue  # Skip the header

                fields = line.strip().split('\t')

                # Check that the necessary fields are numeric before proceeding
                try:
                    start = int(fields[1])  # start_pos
                    end = int(fields[2])    # end_pos
                except ValueError:
                    print(f"Skipping invalid line {i}: {line.strip()}")
                    continue  # Skip lines where start or end position is not numeric

                chrom = fields[0]          # chr_id
                entrez = fields[3]         # entrez_gene_id
                prot_change = fields[4]    # prot_change

                # Find matching annotations
                if chrom in annotations:
                    for feature_start, feature_end, annotation in annotations[chrom]:
                        if start >= feature_start and end <= feature_end:
                            ensp = None
                            for field in annotation.split(';'):
                                if 'protein_id=' in field:
                                    ensp = field.split('=')[1]
                                    break
                            if ensp:
                                # Add match to batch with renamed fields
                                batch.append([chrom, start, end, entrez, prot_change, ensp])

                # Write batch to file every 'batch_size' records
                if len(batch) >= batch_size:
                    writer.writerows(batch)
                    batch.clear()  # Clear batch after writing to file
                    print(f"Processed {i} lines so far...")  # Status update

            # Write remaining entries in the batch
            if batch:
                writer.writerows(batch)
                print("Wrote remaining records to file.")

    print(f"Process completed. Results written to {output_file}")


# Run the workflow
parse_and_serialize_ann()  # Run once to create the serialized annotations file if needed
process_large_input_bed()  # Process large 'input_bed' and write results
