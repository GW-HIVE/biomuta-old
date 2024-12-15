import csv
import pandas as pd
import json
import logging
import os
import glob
import time

# Logging
logging.basicConfig(filename="combine_cbio2.log",
                    filemode='a',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logging.info("Logger started ----------------------")

# Paths
base_csv_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/chr_pos_to_ensp.tsv'
json_dir_path = '/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations_toy'
study_ids_with_do_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/study_ids_with_do.json'

# Debugging
import traceback

'''
# Check for unusually large fields
# Check for rows with more than 6 columns in base_csv_path
with open(base_csv_path, 'r') as file:
    reader = csv.reader(file, delimiter="\t")
    i = 0
    try:
        for i, row in enumerate(reader):
            if len(row) != 6: # Expected number of columns is 6
                logging.warning(f"Row {i + 1}: {row}") # Print problematic row
            for field in row:
                if len(field) > 131072: # Size larger than the default limit
                    logging.warning(f"Large field at row {i + 1}: {field[:100]}...") # Print first 100 chars for preview
    except Exception as e:
        logging.warning(f"Problematic index: {i}")
        logging.error(e)
        logging.error(traceback.format_exc)

'''
### Load base CSV into a dictionary for quick lookup
start_time = time.time()
logging.info("Starting to read the base CSV in chunks...")

# Initialize base_dict
base_dict = {}

# Define chunk size
chunk_size = 1000  # Adjust based on system's memory and processing speed

# Read the CSV file in chunks
for chunk_index, chunk in enumerate(
    pd.read_csv(
        base_csv_path,
        sep='\t',
        dtype={'chr_id': str, 'entrez_gene_id': str, 'prot_change': str},
        low_memory=False,
        chunksize=chunk_size,
    )
):
    logging.info(f"Processing chunk {chunk_index + 1}...")

    # Add chunk data to the dictionary
    for _, row in chunk.iterrows():
        key = (str(row['chr_id']), str(row['entrez_gene_id']), str(row['prot_change']))
        base_dict[key] = row

    logging.info(f"Chunk {chunk_index + 1} processed. Current dict size: {len(base_dict)}")

logging.info(f"Finished processing CSV. Total keys in base_dict: {len(base_dict)}")
logging.info(f"Base dictionary created in {time.time() - start_time:.2f} seconds.")





''' This takes too long
base_df = pd.read_csv(
    base_csv_path, sep='\t', dtype={'chr_id': str, 'entrez_gene_id': str, 'prot_change': str}, low_memory=False
)
logging.info(f"Base CSV loaded with {len(base_df)} rows.")

# Build the dictionary with updated chromosome IDs
base_dict = {
    (str(row['chr_id']), str(row['entrez_gene_id']), str(row['prot_change'])): row # Convert all to string to avoid unexpected dtypes
    for _, row in base_df.iterrows()
}
'''

# Check the contents of the first 10 pairs and their dtypes
for i, (key, value) in enumerate(base_dict.items()):
    if i >= 10:
        break
    logging.info(f"Key: {key} | Key Types: {tuple(type(k) for k in key)}")
    logging.info(f"Value: {value.to_dict()}") # Print value as a dictionary for readability
    value_types = {col: type(value[col]) for col in value.index}
    logging.info(f"Value Types: {value_types}")
    logging.info("-" * 50)

# Load study IDs with DO information
with open(study_ids_with_do_path) as f:
    study_ids_with_do = {entry["studyId"]: entry["do_name"] for entry in json.load(f)}

# Collect output data in a list
output_data =[]

# Loop through JSON files in the mutation directory
json_files = glob.glob(os.path.join(json_dir_path, '*.json'))
total_files = len(json_files)
processed_mutations = 0

# Debugging json processing
processed_rows = []

for i, json_file in enumerate(json_files, start=1):
    with open(json_file) as f:
        mutations = json.load(f)
        logging.info(f"Processing file {i}/{total_files}: {json_file}")
        for mutation in mutations:
            key = (str(mutation['chr']), str(mutation['entrezGeneId']), str(mutation['proteinChange'])) #check what is stored
            #dict is a hash map in python. How does key hashing work for tuples or -> .join() to one string
            # Retrieve the corresponding row from base_dict
            row = base_dict.get(key) #print row

            # Store the processed row for debugging
            if len(processed_rows) < 10:
                processed_rows.append((key, row))

# Print they first 10 rows and their data types
for i, (key, row) in enumerate(processed_rows):
    logging.info(f"Row {i + 1}:")
    logging.info(f"Key: {key} | Key Types: {tuple(type(k) for k in key)}")
    if row is None:
        logging.warning(f"Row is None, skipping...")
        logging.info("-" * 50)
        continue
    logging.info(f"Row data: {row.to_dict()}") # Convert row (Pandas Series) to dict for readability
    value_types = {col: type(row[col]) for col in row.index}
    logging.info(f"Row value types: {value_types}")
    logging.info("-" * 50)

'''
            if row:
                # Extract data for the new CSV
                sample_name = mutation.get('sampleId', '')
                ref_nt = mutation.get('referenceAllele', '')
                alt_nt = mutation.get('variantAllele', '')
                aa_pos = mutation.get('proteinPosStart', '')
                ref_aa = row['prot_change'][0]
                alt_aa = row['prot_change'][-1]
                do_name = study_ids_with_do.get(mutation['studyId'], '')
                uniprotkb_canonical_ac = row['uniprot_canonical_ac']

                # Append data to the output list
                output_data.append({
                    "sample_name": sample_name,
                    "chr_id": row['chr_id'],
                    "start_pos": row['start_pos'],
                    "end_pos": row['end_pos'],
                    "ref_nt": ref_nt,
                    "alt_nt": alt_nt,
                    "aa_pos": aa_pos,
                    "ref_aa": ref_aa,
                    "alt_aa": alt_aa,
                    "do_name": do_name,
                    "uniprotkb_canonical_ac": uniprotkb_canonical_ac,
                    "source": "cbioportal",
                    "dbsnp_id": ''
                })

            # Increment mutation counter
            processed_mutations += 1

            # Print progress every 1000 mutations
            if processed_mutations % 1000 == 0:
                logging.info(f"Processed {processed_mutations} mutations...")

    # Print progress after each file is processed
    logging.info(f"Completed {i}/{total_files} files ({json_file})")

# Convert collected data to DataFrame and save to CSV
output_df = pd.DataFrame(output_data)
output_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/final_output.csv'
output_df.to_csv(output_path, index=False)

logging.info(f"CSV file saved to {output_path}")
'''