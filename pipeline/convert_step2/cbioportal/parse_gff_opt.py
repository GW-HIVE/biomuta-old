'''
This is the optimized version of parse_gff.py. It hasn't been tested yet. The non-optimized parse_gff.py took over 48h to run.
'''
import csv
import gffutils
import logging
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
from utils import ROOT_DIR
from utils.config import get_config

logging.basicConfig(filename="cancer_mapping.log",
                    filemode='a',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logging.info("Logger started ----------------------")

# Create a database from gff3 (on-time operation)
'''
This will parse the file, infer the relationships among
the features in the file, and store the features and
relationships in the database file.
'''
#gff_file = '/data/shared/repos/biomuta-old/downloads/ensembl/Homo_sapiens.GRCh38.113.gff3'
database = '/data/shared/repos/biomuta-old/downloads/ensembl/Homo_sapiens.GRCh38.113.db'
'''
db = gffutils.create_db(
    gff_file,
    database,
    id_spec=['ID', 'Name', 'Parent'],   # Try using ID, then Name, then Parent
    merge_strategy='create_unique',     # Append unique suffixes to duplicate IDs
    force=True,                         # Overwrite if needed
    keep_order=True                     # Maintain original order
    )
'''

# Load the database
db = gffutils.FeatureDB(database)

'''
# Get all unique feature types
feature_types = set()
for feature in db.all_features():
    feature_types.add(feature.featuretype)

logging.info(f"Available feature types: {feature_types}")
'''

def get_ensp_for_position(chrom, start, end):
    '''
    Query the GFF database to find ENSP IDs for the given chromosomal position.
    Args:
        chrom (str): Chromosome number, e.g. 'chr1' or 'chrX'
        start (int): Start position
        end (int): End position
    Returns:
        list: List of ENSP IDs (or empty list if none found)
    '''
    ensp_ids = []

    # Adjust chromosome format if needed
    if chrom.startswith('chr'):
        chrom = chrom.replace('chr', '')

    # Adjust start to 1-based
    start += 1

    # Query the database for 'CDS' features in the given region
    for feature in db.region(region=(chrom, start, end), featuretype='CDS'):
        if 'protein_id' in feature.attributes:
            ensp_ids.extend(feature.attributes['protein_id'])
    return ensp_ids if ensp_ids else ['N/A']


# Helper function to clean and map chromosome IDs back to their original format
chr_id_cache = {} # Store IDs in a dictionary to avoid repeated computations, esp helpful since input has many repeated chr IDs
def clean_chr_id(chr_id):
    if chr_id not in chr_id_cache:
        chr_id_clean = chr_id.lstrip("chr") # Remove 'chr' prefix
        chr_id_clean = "23" if chr_id_clean == "X" else "24" if chr_id_clean == "Y" else chr_id_clean
        chr_id_cache[chr_id] = chr_id_clean
    return chr_id_cache[chr_id]


def process_bed_file(input_file, output_file):
    '''
    Process a BED file (TSV format) to get ENSP IDs for each position and save the results.
    '''
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Write header to output file
        writer.writerow(['chr_id', 'start_pos', 'end_pos', 'entrez_gene_id', 'prot_change', 'ensp'])

        for row in reader:
            if not row or row[0].startswith('chr_id'):
                continue # Skip header or empty rows

            chr_id, start_pos, end_pos, entrez_gene_id, prot_change = row
            start_pos, end_pos = int(start_pos), int(end_pos) # Ensure positions are integers

            ensp_ids = get_ensp_for_position(chr_id, start_pos, end_pos) # Fetch ENSP IDs for the given position

            # Write each ENSP as a separate row if there are more than 1
            if ensp_ids:
                for ensp in ensp_ids:
                    writer.writerow([clean_chr_id(chr_id), start_pos, end_pos, entrez_gene_id, prot_change, ensp])

# Example usage
'''
chrom = '10'
start = 43163970
end = 43163971
ensps = get_ensp_for_position(chrom, start, end)
logging.info(f"ENSP IDs: {ensps}")
'''
config_obj = get_config()
wd = Path(config_obj["relevant_paths"]["generated_datasets"])
input_file = wd / '2024_10_22' / 'liftover' / 'hg38_combined.bed'  # Write a util to get latest dir
output_file = wd / '2024_10_22' / 'mapping_ids' / 'chr_pos_to_ensp.tsv'

process_bed_file(input_file, output_file)