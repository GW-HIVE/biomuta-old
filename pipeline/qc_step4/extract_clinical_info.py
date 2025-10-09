import glob
import json
import logging
import os
import pandas as pd

# Logging
logging.basicConfig(filename="extract_clinical_info.log",
                    filemode='a',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logging.info("Logger started ----------------------")

# Paths
base_dict_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/base_dict.json'
json_dir_path = '/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations'
bm_path = "/data/shared/repos/biomuta-old/generated_datasets/current/final_cbio_table.csv"
output_path = '/data/shared/repos/biomuta-old/generated_datasets/current/clinical-information.csv'

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

df = pd.read_csv(bm_path)

# Collect output data in a list
output_data = []

# Loop through JSON files in the mutation directory
json_files = glob.glob(os.path.join(json_dir_path, '*.json'))
total_files = len(json_files)
processed_mutations = 0
batch_size = 1000000  # Number of mutations to process before writing to CSV
batch_counter = 0

# Define a temporary output directory for intermediate CSVs
temp_dir = "/data/shared/repos/biomuta-old/generated_datasets/temp_clinical_outputs"
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
                sample_id = mutation.get('sampleId', '')
                patient_id = mutation.get('patientId', '')
                study_id = mutation.get('studyId', '')
                ref_nt = mutation.get('referenceAllele', '')
                chr_id = mutation.get('chr', '')
                alt_nt = mutation.get('variantAllele', '')
                aa_pos = mutation.get('proteinPosStart', '')

                # Append data to the output list
                output_data.append({
                    "sample_id": sample_id,
                    "patient_id": patient_id,
                    "proj__project_id": study_id,
                    "ref_nt": ref_nt,
                    "chr_id": chr_id,
                    "alt_nt": alt_nt,
                    "aa_pos": aa_pos
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