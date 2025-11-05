import os
import json

# Directory containing JSON files
json_files = "/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations"

# Iterate over each file in the directory
for filename in os.listdir(json_files):
    if filename.endswith(".json"):  # Make sure it's a JSON file
        file_path = os.path.join(json_files, filename)
        
        with open(file_path, 'r') as f:
            data = json.load(f)
            
            # Check if "proteinChange" key exists and its value is an empty string
            if "proteinChange" in data and data["proteinChange"] == "":
                print(f"File {filename} contains an empty 'proteinChange' value.")
                break  # Stop at the first occurrence
