import pandas as pd
import json
import os

# Initialize the lookup dictionary
mutation_lookup = {}

# Path to JSON files
json_directory = '/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations'

# Build the lookup table from JSON files
for json_file in os.listdir(json_directory):
    if json_file.endswith('.json'):
        with open(os.path.join(json_directory, json_file), 'r') as f:
            json_data = json.load(f)
        for entry in json_data:
            chr_id = entry.get('chr')
            key = (chr_id, entry.get('entrezGeneId'), entry.get('proteinChange'))
            mutation_lookup[key] = entry.get('studyId')

# Process the CSV file in chunks
csv_file_path = '/data/shared/biomuta/generated/datasets/2024_10_22/mapping_ids/chr_pos_to_ensp_old.csv'
updated_csv_file_path = '/data/shared/biomuta/generated/datasets/2024_10_22/mapping_ids/chr_pos_to_ensp_old_study_id.csv'

chunksize = 10000  # Adjust chunk size based on memory availability
output_chunks = []

# Function to handle matching with debugging and special case handling for X and Y
def debug_matching(row):
    # Skip rows where 'prot_change' is empty or NaN
    if pd.isna(row['prot_change']) or row['prot_change'] == '':
        return ''  # Skip the row and don't add a studyId
    
    # Remove "chr" prefix from chr_id for matching
    chr_stripped = row['chr_id'].replace('chr', '')
    
    # Handle special cases for chrX and chrY
    if chr_stripped == 'X':
        chr_stripped = '23'
    elif chr_stripped == 'Y':
        chr_stripped = '24'
    
    key = (chr_stripped, row['entrez'], row['prot_change'])
    
    if key in mutation_lookup:
        return mutation_lookup[key]
    else:
        print(f"No match for key: {key}")  # Debugging output
        return ''

# Process the CSV file in chunks with the updated key logic
output_chunks = []
for chunk in pd.read_csv(csv_file_path, chunksize=chunksize):
    # Apply the matching function and skip rows where 'prot_change' is empty
    chunk['studyId'] = chunk.apply(debug_matching, axis=1)
    output_chunks.append(chunk)

# Save the final DataFrame to a new CSV file
final_df = pd.concat(output_chunks)
final_df.to_csv(updated_csv_file_path, index=False)
