import csv
import gffutils
import logging
import sys
from pathlib import Path
from collections import defaultdict
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
from utils import ROOT_DIR
from utils.config import get_config

logging.basicConfig(filename="/data/shared/repos/biomuta-old/pipeline/logs/2_convert/cancer_mapping.log",
                    filemode='a',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logging.info("Logger started ----------------------")

# Load the database
database = '/data/shared/repos/biomuta-old/downloads/ensembl/Homo_sapiens.GRCh38.113.db'
db = gffutils.FeatureDB(database)

def build_cds_lookup_table():
    """
    Pre-build a lookup table of all CDS features. 
    This runs once and creates an in-memory structure for fast lookups.
    Expected time: 5-15 minutes for the entire genome.
    """
    logging.info("Building CDS lookup table - this may take 10-15 minutes but will save hours later...")
    
    # Structure: {chromosome: [(start, end, protein_ids), ...]}
    cds_lookup = defaultdict(list)
    
    feature_count = 0
    for feature in db.features_of_type('CDS'):
        if 'protein_id' in feature.attributes:
            chrom = str(feature.seqid)
            # Remove 'chr' prefix if present for consistency
            if chrom.startswith('chr'):
                chrom = chrom[3:]
                
            protein_ids = feature.attributes['protein_id']
            cds_lookup[chrom].append((feature.start, feature.end, protein_ids))
            
        feature_count += 1
        if feature_count % 50000 == 0:
            logging.info(f"Processed {feature_count} CDS features...")
    
    # Sort intervals by start position for each chromosome for faster searching
    for chrom in cds_lookup:
        cds_lookup[chrom].sort(key=lambda x: x[0])
        logging.info(f"Chromosome {chrom}: {len(cds_lookup[chrom])} CDS features")
    
    logging.info(f"CDS lookup table built with {feature_count} total features")
    return dict(cds_lookup)

def find_overlapping_cds(cds_lookup, chrom, start, end):
    """
    Fast lookup using sorted intervals.
    Much faster than database queries.
    """
    # Clean chromosome format
    if chrom.startswith('chr'):
        chrom = chrom[3:]
    
    if chrom not in cds_lookup:
        return ['N/A']
    
    # Convert to 1-based coordinates to match GFF
    start += 1
    
    protein_ids = []
    
    # Binary search could be implemented here, but linear search should be fast enough
    # for most chromosomes given the sorting
    for cds_start, cds_end, prot_ids in cds_lookup[chrom]:
        # Skip if CDS is completely before our interval
        if cds_end < start:
            continue
        # Stop if CDS is completely after our interval (since list is sorted)
        if cds_start > end:
            break
        # Check for overlap
        if cds_start <= end and cds_end >= start:
            protein_ids.extend(prot_ids)
    
    return protein_ids if protein_ids else ['N/A']

def clean_chr_id(chr_id):
    """Clean chromosome ID"""
    chr_id_clean = chr_id.lstrip("chr")
    if chr_id_clean == "X":
        return "23"
    elif chr_id_clean == "Y":
        return "24"
    return chr_id_clean

def process_bed_file_fast(input_file, output_file, cds_lookup):
    """
    Process BED file using pre-built lookup table.
    Expected time: 30-60 minutes for 5.4M rows.
    """
    logging.info("Starting fast BED file processing...")
    
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        # Write header
        writer.writerow(['chr_id', 'start_pos', 'end_pos', 'entrez_gene_id', 'prot_change', 'ensp'])
        
        processed = 0
        for row in reader:
            if not row or row[0].startswith('chr_id'):
                continue
            
            chr_id, start_pos, end_pos, entrez_gene_id, prot_change = row
            start_pos, end_pos = int(start_pos), int(end_pos)
            
            # Fast lookup using pre-built table
            protein_ids = find_overlapping_cds(cds_lookup, chr_id, start_pos, end_pos)
            
            # Write results
            for protein_id in protein_ids:
                writer.writerow([clean_chr_id(chr_id), start_pos, end_pos, entrez_gene_id, prot_change, protein_id])
            
            processed += 1
            if processed % 100000 == 0:
                logging.info(f"Processed {processed:,} rows...")
    
    logging.info(f"Completed processing {processed:,} rows")

def process_bed_file_pandas(input_file, output_file, cds_lookup):
    """
    Alternative pandas-based approach for even better performance.
    Expected time: 15-30 minutes for 5.4M rows.
    """
    logging.info("Loading BED file with pandas...")
    
    # Read the BED file
    df = pd.read_csv(input_file, sep='\t', header=0, 
                     names=['chr_id', 'start_pos', 'end_pos', 'entrez_gene_id', 'prot_change'])
    
    logging.info(f"Loaded {len(df):,} rows")
    
    # Process in chunks to manage memory
    chunk_size = 50000
    results = []
    
    for i in range(0, len(df), chunk_size):
        chunk = df.iloc[i:i+chunk_size]
        logging.info(f"Processing chunk {i//chunk_size + 1}/{(len(df)//chunk_size) + 1}")
        
        chunk_results = []
        for _, row in chunk.iterrows():
            protein_ids = find_overlapping_cds(cds_lookup, row['chr_id'], 
                                             row['start_pos'], row['end_pos'])
            
            for protein_id in protein_ids:
                chunk_results.append({
                    'chr_id': clean_chr_id(row['chr_id']),
                    'start_pos': row['start_pos'],
                    'end_pos': row['end_pos'],
                    'entrez_gene_id': row['entrez_gene_id'],
                    'prot_change': row['prot_change'],
                    'ensp': protein_id
                })
        
        results.extend(chunk_results)
    
    # Save results
    result_df = pd.DataFrame(results)
    result_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Saved {len(result_df):,} results to {output_file}")

# Main execution
if __name__ == "__main__":
    config_obj = get_config()
    wd = Path(config_obj["relevant_paths"]["generated_datasets"])
    input_file = wd / 'current' / 'liftover' / 'hg38_combined.bed'
    output_file = wd / 'current' / 'mapping_ids' / 'chr_pos_to_ensp.tsv'
    
    # Create output directory if it doesn't exist
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Build lookup table (one-time cost: ~10-15 minutes)
    cds_lookup = build_cds_lookup_table()
    
    # Step 2: Process BED file using lookup table (~30-60 minutes)
    # Choose one approach:
    
    # Option A: Simple CSV approach
    process_bed_file_fast(input_file, output_file, cds_lookup)
    
    # Option B: Pandas approach (potentially faster)
    # process_bed_file_pandas(input_file, output_file, cds_lookup)
    
    logging.info("Script completed successfully!")
    logging.info("Script finished ----------------------")