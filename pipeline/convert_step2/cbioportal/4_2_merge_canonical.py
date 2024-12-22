# Concatenated this script to 4_1_canonical_yes_no.py

import json

# File paths
file1 = "ensp_to_uniprot_canonical.json"
file2 = "ensp_to_uniprot_canonical_gffutils.json"
output_file = "merged_ensp_to_uniprot.json"

# Load JSON data
with open(file1, "r") as f:
    data1 = json.load(f)

with open(file2, "r") as f:
    data2 = json.load(f)

# Filter out entries where "canonical" == "no"
filtered_data1 = {k: v for k, v in data1.items() if v.get("canonical") == "yes"}
filtered_data2 = {k: v for k, v in data2.items() if v.get("canonical") == "yes"}

# Merge dictionaries
merged_data = {**filtered_data1, **filtered_data2}

# Remove the "canonical" key from all entries
for entry in merged_data.values():
    entry.pop("canonical", None)

# Flatten the JSON structure
flattened_data = {key: value['uniprotkb_canonical_ac'] for key, value in merged_data.items()}

# Save the cleaned, merged data to a new JSON file
with open(output_file, "w") as f:
    json.dump(flattened_data, f, indent=4)

print(f"Merged data saved to {output_file}")
