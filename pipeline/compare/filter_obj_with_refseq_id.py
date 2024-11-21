import json

# List of refseqMrnaId substrings to match
refseq_ids_to_filter = ["NM_005228", "NM_201282", "NM_201283", "NM_201284", "NM_001346897", "NM_001346898", "NM_001346899", "NM_001346900", "NM_001346941"]

# Read the list of JSON files from 'files_with_refseq_ids.txt'
with open('/data/shared/biomuta/generated/datasets/2024_10_22/files_with_refseq_ids.txt', 'r') as file_list:
    json_files = file_list.readlines()

# Initialize a list to store the filtered JSON objects
filtered_data = []

# Loop through each file listed in 'files_with_refseq_ids.txt'
for json_file in json_files:
    json_file = json_file.strip()  # Remove any leading/trailing whitespaces or newlines
    file_path = f'/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations/{json_file}'
    
    # Open and read the JSON data from each file
    with open(file_path, 'r') as f:
        data = json.load(f)
        
        # Filter the JSON objects based on refseqMrnaId
        for record in data:
            if any(refseq_id in record['refseqMrnaId'] for refseq_id in refseq_ids_to_filter):
                filtered_data.append(record)

# Write the filtered data to a new JSON file
output_file = '/data/shared/biomuta/generated/datasets/2024_10_22/filtered_by_refseq_id.json'
with open(output_file, 'w') as f:
    json.dump(filtered_data, f, indent=4)

print(f"Filtered data has been written to {output_file}")
