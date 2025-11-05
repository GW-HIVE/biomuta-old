import pandas as pd
import json

# Initialize the lookup dictionary
do_lookup = {}

# Path to the DO mapping JSON file
do_map = '/data/shared/biomuta/generated/datasets/2024_10_22/study_ids_with_do.json'

# Build the lookup dictionary
with open(do_map, 'r') as f:
    study_id_do = json.load(f)
    for entry in study_id_do:
        key = entry.get('studyId')
        do_lookup[key] = entry.get('do_name')

# Process the CSV file
csv_file_path = '/data/shared/biomuta/generated/datasets/2024_10_22/mapping_ids/chr_pos_to_ensp_old_study_id.csv'
updated_csv_file_path = '/data/shared/biomuta/generated/datasets/2024_10_22/mapping_ids/chr_pos_to_ensp_old_do_name.csv'

# Load the CSV into a DataFrame
df = pd.read_csv(csv_file_path)

# Add a new column 'do_name' by mapping 'studyId' using the lookup dictionary
df['do_name'] = df['studyId'].apply(lambda x: do_lookup.get(x, ''))

# Save the updated DataFrame to a new CSV file
df.to_csv(updated_csv_file_path, index=False)

print(f"Output CSV with 'do_name' written to {updated_csv_file_path}")
