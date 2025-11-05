import json

# Path to JSON file
file_path = '/data/shared/biomuta/generated/datasets/2024_10_22/cancer_types_with_do.json'

# Open and load JSON data from the file
with open(file_path, 'r') as file:
    data = json.load(file)

# Extract unique `do_name` values
unique_do_names = set(entry["do_name"] for entry in data)

# Print each `do_name` value on a new line without quotes
for do_name in unique_do_names:
    print(do_name)