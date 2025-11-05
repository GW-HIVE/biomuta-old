import json
import logging
import os
import pandas as pd

# Logging
logging.basicConfig(filename="6_1_create_dict.log",
                    filemode='a',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logging.info("Logger started ----------------------")

# Paths
base_tsv_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/chr_pos_to_ensp.tsv'
uniprot_ac_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/merged_ensp_to_uniprot.json'
base_dict_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/base_dict.json'
chunk_size = 1000

# Initialize base_dict
base_dict = {}

# Load UniProt mapping
with open(uniprot_ac_path, 'r') as f:
    uniprot_mapping = json.load(f)
logging.info(f"UniProt mapping loaded successfully with {len(uniprot_mapping)} entries.")

# Create base_dict
logging.info("Processing base TSV to recreate base_dict...")

# Track logged rows
logged_missing_ensp_keys = set()

# Process chunks
for chunk in pd.read_csv(
    base_tsv_path,
    sep='\t',
    dtype={'chr_id': str, 'entrez_gene_id': str, 'prot_change': str},
    low_memory=False,
    chunksize=chunk_size,
):
    # Create masks for 'prot_change' conditions
    prot_change_mask = chunk['prot_change'].apply(lambda x: x == "MUTATED" or x.startswith("*"))
    # Identify rows with missing 'ensp'
    missing_ensp_mask = chunk['ensp'].isna() | (chunk['ensp'] == "") | (chunk['ensp'].apply(lambda x: not x))
    # Filter rows with missing 'ensp' and valid 'prot_change'
    missing_ensp_valid_rows = chunk[missing_ensp_mask & ~prot_change_mask]
    # Log missing 'ensp' for valid rows only (those that do not have prot_change == "MUTATED" or starting with "*")
    for _, row in missing_ensp_valid_rows.iterrows():
        # Create the key
        keys = ["chr_id", "entrez_gene_id", "prot_change"]
        list_key = [row[x] for x in keys]
        key = "|".join(list_key)
        if key not in logged_missing_ensp_keys:
            logging.warning(f"Row with key {key} has 'ensp' as None or empty (NaN)")
            logging.warning(f"Row: {row}")
            logged_missing_ensp_keys.add(key) # Mark this key as logged

    # Filter out rows with missing 'ensp' for further processing
    valid_rows = chunk[~missing_ensp_mask]
    # Process valid rows
    for _, row in chunk.iterrows(): # type(row) == pd.Series
        # Skip rows where prot_change is invalid
        if prot_change_mask.get(_, False):
            continue # Skip these rows
        # Create the key
        keys = ["chr_id", "entrez_gene_id", "prot_change"]
        list_key = [row[x] for x in keys]
        key = "|".join(list_key)
        # Check if 'ensp' is missing
        ensp_id = str(row.get('ensp', "")).strip()
        if not ensp_id:
            continue # Skip this row if 'ensp' is missing
        
        # Map ENSP to UniProt
        uniprot_id = uniprot_mapping.get(ensp_id, None)
        if uniprot_id is None:
            continue  # Skip this row

        # Check for None or NaN in uniprot_mapping keys
        for ensp_id, uniprot_id in uniprot_mapping.items():
            if pd.isna(uniprot_id):
                logging.warning(f"Found NaN or None for ENSP ID {ensp_id} in UniProt mapping.")

        # Replace 'ensp' with 'uniprot_canonical_ac'
        if uniprot_id:
            row['uniprot_canonical_ac'] = uniprot_id
        else:
            row['uniprot_canonical_ac'] = None
        # Drop the original 'ensp' column
        row = row.drop(labels=['ensp'])

        # Add the row as a value in the base dictionary
        base_dict[key] = row.to_list()

logging.info(f"base_dict contains {len(base_dict)} keys before saving.")

# Save base_dict to JSON
with open(base_dict_path, "w") as f:
    json.dump(base_dict, f, indent=4)
logging.info(f"base_dict successfully saved to {base_dict_path}")

for i, (key, value) in enumerate(base_dict.items()):
    logging.debug(f"Key {i}: {key}, Value: {value}")
    if i >= 5:  # Log only the first 5 entries
        break
if len(base_dict) == 0:
    logging.error("base_dict is empty. Skipping save operation.")
else:
    try:
        with open(base_dict_path, "w") as f:
            json.dump(base_dict, f, indent = 4)
        logging.info("base_dict successfully saved.")
    except Exception as e:
        logging.error(f"Failed to save base_dict: {e}")

# Verify the file was saved correctly
if os.path.getsize(base_dict_path) > 0:
    logging.info(f"base_dict successfully saved to {base_dict_path}.")
else:
    logging.error(f"Failed to save base_dict to {base_dict_path}. File is empty.")

# Later when you load in your json file
# Here starts the for loop
#row = base_dict.get(key)
#if row is not None:
#    ref_aa = row[2][0] # index of prot_change is e.g. 2, and the first character is index 0
# ... rest of the fields

# Check the contents of the first 10 pairs and their dtypes
for i, (key, value) in enumerate(base_dict.items()):
    if i >= 10:
        break
    logging.info(f"Key: {key} | Key Type: {type(key)}")
    logging.info(f"Value: {value} | Value Type: {type(value)}")
    logging.info("-" * 50)