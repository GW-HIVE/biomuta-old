# Should I keep rows whose only difference with each other is sample name?

# QC
#   - How many mutations did I lose at each step in the pipeline, e.g.
#      how many weren't lifted over
#      how many were isoforms (ENSP IDs thrown away because they are isoforms)
#      how many weren't found in UniProt map

# Don't strip the number after uniprot canonical ac (upstream in the pipeline)
#   1. Re-run 3_1_ensp_to_uniprot.py: Done! (But there was no need because I updated base_dict directly. Useful for the next release.)
#   2. Re-run 3_2_ensp_to_uniprot.sh: No need because the API doesn't specify the canonical isoform (e.g. P29692). Could map the resulting accessions to masterlist to extract (e.g. P29692-1). Updating base_dict... Done!
#   3. Use base_dict_updated to rebuild the final table. Add quotes and +1 to start_pos: Done!

import glob
import json
import logging
import os
import pandas as pd

# Logging
logging.basicConfig(filename="/data/shared/repos/biomuta-old/pipeline/logs/3_combine/6_combine_cbio6.log",
                    filemode='a',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logging.info("Logger started ----------------------")

# Paths
base_dict_path = '/data/shared/repos/biomuta-old/generated_datasets/current/base_dict.json'
json_dir_path = '/data/shared/biomuta/downloads/cbioportal/current/mutations'
study_ids_with_do_path = '/data/shared/repos/biomuta-old/generated_datasets/current/study_ids_with_do.json'
output_path = '/data/shared/repos/biomuta-old/generated_datasets/current/final_cbio_table.csv'

# Load base_dict from file
if os.path.exists(base_dict_path):
    if os.path.getsize(base_dict_path) > 0:
        logging.info("base_dict found. Attempting to load it...")
        with open(base_dict_path, 'r') as f:
            base_dict = json.load(f)
        logging.info("base_dict loaded successfully.")
    else:
        logging.error(f"The file {base_dict_path} is empty.")
else:
    logging.error(f"The file {base_dict_path} does not exist.")

# Load study IDs with DO information
with open(study_ids_with_do_path) as f:
    study_ids_with_do = {entry["studyId"]: entry["do_name"] for entry in json.load(f)}

# Collect output data in a list
output_data = []

# Loop through JSON files in the mutation directory
json_files = glob.glob(os.path.join(json_dir_path, '*.json'))
total_files = len(json_files)
processed_mutations = 0
batch_size = 1000000  # Number of mutations to process before writing to CSV
batch_counter = 0

# Define a temporary output directory for intermediate CSVs
temp_dir = "/data/shared/repos/biomuta-old/generated_datasets/temp_outputs"
os.makedirs(temp_dir, exist_ok=True)

for i, json_file in enumerate(json_files, start=1):
    with open(json_file) as f:
        mutations = json.load(f)
        logging.info(f"Processing file {i}/{total_files}: {json_file}")
        for mutation in mutations:
            keys = (str(mutation['chr']), str(mutation['entrezGeneId']), str(mutation['proteinChange']))
            logging.debug(f"Keys: {keys}")
            list_key = list(keys)
            key = "|".join(list_key)
            logging.debug(f"Joined key: {key}")
            # Retrieve the corresponding row from base_dict
            row = base_dict.get(key)
            logging.debug(f"Row structure: {row}")

            if row is not None: # Check if a row exists for the key
                # Extract data for the new CSV
                sample_name = mutation.get('sampleId', '')
                ref_nt = mutation.get('referenceAllele', '')
                alt_nt = mutation.get('variantAllele', '')
                aa_pos = mutation.get('proteinPosStart', '')
                ref_aa = row[4][0] # index of prot_change is e.g. 2, and the first character is index 0
                alt_aa = row[4][-1]
                do_name = study_ids_with_do.get(mutation['studyId'], '')
                uniprotkb_canonical_ac = row[5]

                # Append data to the output list
                output_data.append({
                    "sample_name": sample_name,
                    "chr_id": row[0],
                    "start_pos": row[1] + 1,
                    "end_pos": row[2],
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

            # Write output_data to CSV in batches
            if len(output_data) >= batch_size:
                batch_counter += 1
                batch_path = os.path.join(temp_dir, f"batch_{batch_counter}.csv")
                pd.DataFrame(output_data).to_csv(batch_path, index=False)
                logging.info(f"Batch {batch_counter} written to {batch_path}")
                output_data = []  # Clear the list to free memory

        logging.info(f"Completed {i}/{total_files} files ({json_file})")

# Write remaining data after the loop
if output_data:
    batch_counter += 1
    batch_path = os.path.join(temp_dir, f"batch_{batch_counter}.csv")
    pd.DataFrame(output_data).to_csv(batch_path, index=False)
    logging.info(f"Final batch {batch_counter} written to {batch_path}")

# Combine all temporary CSVs into the final output
temp_files = glob.glob(os.path.join(temp_dir, "*.csv"))
combined_df = pd.concat([pd.read_csv(file) for file in temp_files], ignore_index=True)
combined_df.to_csv(output_path, index=False, quoting=1)
logging.info(f"Final combined CSV saved to {output_path}")
