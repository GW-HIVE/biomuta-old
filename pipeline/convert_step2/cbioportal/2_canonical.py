import json
import pandas as pd

# Load the ENSP to UniProt mapping JSON
with open("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/ensp_to_uniprot_mappings_toy.json", "r") as f:
    ensp_to_uniprot = json.load(f)

# Load the isoform data CSV
isoform_data = pd.read_csv("/data/shared/repos/biomuta-old/downloads/glygen/human_protein_masterlist.csv")

# Prepare a dictionary to store the results
result = {}

# Function to strip suffixes (anything after a hyphen) from both isoform IDs and UniProtKB Canonical ACs
def strip_suffix(identifier):
    if isinstance(identifier, str) and '-' in identifier:
        return identifier.split('-')[0]  # Strip everything after the first hyphen
    return identifier

# Iterate over each ENSP and its corresponding UniProt ID
for ensp, uniprot in ensp_to_uniprot.items():
    # Check for matching UniProt IDs in either reviewed_isoforms or unreviewed_isoforms
    for _, entry in isoform_data.iterrows():
        # Strip suffixes from isoform IDs before comparison
        reviewed_isoforms = strip_suffix(entry.get("reviewed_isoforms", ""))
        unreviewed_isoforms = strip_suffix(entry.get("unreviewed_isoforms", ""))
        
        # Check if the UniProt ID matches any isoform (stripped version)
        if uniprot == reviewed_isoforms or uniprot == unreviewed_isoforms:
            # If a match is found, add the uniprotkb_canonical_ac to the result
            uniprotkb_ac = strip_suffix(entry.get("uniprotkb_canonical_ac"))
            # Store the first match found for each ENSP identifier
            result[ensp] = uniprotkb_ac
            break  # Exit inner loop once the first match is found

# Write the result to a JSON file
with open("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/canonical_toy.json", "w") as json_file:
    json.dump(result, json_file, indent=4)
