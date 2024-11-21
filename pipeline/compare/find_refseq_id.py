import os
import json

# Set the directory containing your JSON files and the list of refseqMrnaIds to search for
directory = '/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations'
refseq_ids_to_find = ["NM_005228", "NM_201282", "NM_201283", "NM_201284", "NM_001346897", "NM_001346898", "NM_001346899", "NM_001346900", "NM_001346941"]  # Add more IDs as needed
output_file = '/data/shared/biomuta/generated/datasets/2024_10_22/files_with_refseq_ids.txt'

# Open the output file in write mode
with open(output_file, 'w') as outfile:
    # Loop through all files in the specified directory
    for filename in os.listdir(directory):
        if filename.endswith('.json'):
            file_path = os.path.join(directory, filename)
            # Open and load each JSON file
            with open(file_path, 'r') as json_file:
                try:
                    data = json.load(json_file)
                    # Check for partial matches with refseqMrnaIds
                    if any(any(refseq_id in item.get("refseqMrnaId", "") for refseq_id in refseq_ids_to_find) for item in data):
                        # Write the filename to the output file if a partial refseqMrnaId is found
                        outfile.write(filename + '\n')
                except json.JSONDecodeError:
                    print(f"Error decoding JSON in file {filename}")
