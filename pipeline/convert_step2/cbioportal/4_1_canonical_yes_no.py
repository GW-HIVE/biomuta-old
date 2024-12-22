# Improvements needed: log IDs that were not found and checksum for the number of input entries vs output

import json
import os
import pandas as pd

# Load the ENSP to UniProt mapping JSON
with open("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/formatted_gffutils_mappings.json", "r") as f:
    ensp_to_uniprot = json.load(f)

# Load the isoform data CSV
isoform_data = pd.read_csv("/data/shared/repos/biomuta-old/downloads/glygen/human_protein_masterlist.csv", usecols=["reviewed_isoforms", "unreviewed_isoforms", "uniprotkb_canonical_ac"])

# Function to strip suffixes (anything after a hyphen) from both isoform IDs and UniProtKB Canonical ACs
def strip_suffix(identifier):
    if isinstance(identifier, str) and '-' in identifier:
        return identifier.split('-')[0]  # Strip everything after the first hyphen
    return identifier

# Preprocess isoform data to remove suffixes
isoform_data["reviewed_isoforms"] = isoform_data["reviewed_isoforms"].apply(strip_suffix)
isoform_data["unreviewed_isoforms"] = isoform_data["unreviewed_isoforms"].apply(strip_suffix)
isoform_data["uniprotkb_canonical_ac"] = isoform_data["uniprotkb_canonical_ac"].apply(strip_suffix)

# Build a lookup dictionary for isoforms
isoform_lookup = {}
for _, row in isoform_data.iterrows():
    for key in ["reviewed_isoforms", "unreviewed_isoforms"]:
        isoform_id = row.get(key, "")
        if isoform_id:
            isoform_lookup[isoform_id] = row.get("uniprotkb_canonical_ac", "")

# Prepare a dictionary to store the results
result = {}

# Iterate over each ENSP and its corresponding UniProt ID
for ensp, uniprot in ensp_to_uniprot.items():
    canonical_ac = isoform_lookup.get(uniprot, None)
    if canonical_ac:
        result[ensp] = {
            "uniprotkb_canonical_ac": canonical_ac,
            "canonical": "yes" if uniprot == canonical_ac else "no"
        }

# Write the result to a JSON file
import os
output_dir = "/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids"
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "ensp_to_uniprot_canonical_gffutils.json")
with open(output_path, "w") as json_file:
    json.dump(result, json_file, indent=4)
