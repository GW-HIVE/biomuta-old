import csv
import pandas as pd
import json
import logging
import os
import glob

# Logging
logging.basicConfig(filename="combine_cbio.log",
                    filemode='a',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

logging.info("Logger started ----------------------")

# Paths
base_csv_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/1-based_chr_pos.csv'
json_dir_path = '/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations_toy'
study_ids_with_do_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/study_ids_with_do.json'

# Debugging

# Load base CSV into a dictionary for quick lookup
# It works
base_df = pd.read_csv(base_csv_path)
base_dict = {
    (str(row['chr_id']), str(row['entrez']), str(row['prot_change'])): row
    for _, row in base_df.iterrows() #might be unexpected data types - make all string
}

'''
# Load study IDs with DO information
with open(study_ids_with_do_path) as f:
    study_ids_with_do = {entry["studyId"]: entry["do_name"] for entry in json.load(f)}

# Collect output data in a list
output_data =[]
'''

# Loop through JSON files in the mutation directory
json_files = glob.glob(os.path.join(json_dir_path, '*.json'))
total_files = len(json_files)
processed_mutations = 0

for i, json_file in enumerate(json_files, start=1):
    with open(json_file) as f:
        mutations = json.load(f)
        for mutation in mutations:
            key = (mutation['chr'], mutation['entrezGeneId'], mutation['proteinChange']) #check what is stored
            logging.info(f"Printing the key...")
            logging.info(key)
            #dict is a hash map in python. How does key hashing work for tuples or -> .join() to one string
            row = base_dict.get(key) #print row
            logging.info(f"Printing the value...")
            logging.info(row)
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