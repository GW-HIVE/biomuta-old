import csv
import json

# File paths
csv_file_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/unique_chr_pos_to_ensp.csv'
json_file_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/canonical.json'
output_file_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/filtered_chr_pos_with_uniprot.csv'

# Load the canonical JSON data
with open(json_file_path, 'r') as json_file:
    canonical_mapping = json.load(json_file)

# Open the CSV file and filter rows based on the canonical JSON keys
with open(csv_file_path, 'r') as csv_file, open(output_file_path, 'w', newline='') as output_file:
    reader = csv.DictReader(csv_file)
    fieldnames = reader.fieldnames + ['uniprot_canonical_ac']  # Add new column
    writer = csv.DictWriter(output_file, fieldnames=fieldnames)
    
    writer.writeheader()  # Write the header to the output file
    
    for row in reader:
        ensp_id = row['ENSP']
        if ensp_id in canonical_mapping:
            # Add the uniprot_canonical_ac value from the JSON data
            row['uniprot_canonical_ac'] = canonical_mapping[ensp_id]
            writer.writerow(row)

print(f"Filtered data written to {output_file_path}")
