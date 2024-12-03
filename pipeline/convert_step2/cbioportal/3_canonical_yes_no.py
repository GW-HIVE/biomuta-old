import json
import pandas as pd
from pathlib import Path

# Load config.json
config_path = Path(__file__).resolve().parent.parent / "config.json"
with open(config_path, "r") as config_file:
    config = json.load(config_file)

# Retrieve paths from updated config
repos_generated_datasets = Path(config["relevant_paths"]["repos_generated_datasets"])
repos_downloads = Path(config["relevant_paths"]["repos_downloads"])
isoform_data_path = repos_downloads / "glygen/human_protein_masterlist.csv"
ensp_to_uniprot_path = repos_generated_datasets / "2024_10_22/mapping_ids/ensp_to_uniprot_mappings_toy.json"
canonical_toy_output_path = repos_generated_datasets / "2024_10_22/mapping_ids/canonical_toy.json"

# Load the ENSP to UniProt mapping JSON
with open(ensp_to_uniprot_path, "r") as f:
    ensp_to_uniprot = json.load(f)

# Load the isoform data CSV
isoform_data = pd.read_csv(isoform_data_path)

# Prepare a dictionary to store the results
result = {}

# Function to strip suffixes (anything after a hyphen) from both isoform IDs and UniProtKB Canonical ACs
def strip_suffix(identifier):
    if isinstance(identifier, str) and '-' in identifier:
        return identifier.split('-')[0]  # Strip everything after the first hyphen
    return identifier

# Iterate over each ENSP and its corresponding UniProt ID
for ensp, uniprot in ensp_to_uniprot.items():
    # Default to "no" for canonical until proven otherwise
    is_canonical = "no"
    
    # Check for matching UniProt IDs in either reviewed_isoforms or unreviewed_isoforms
    for _, entry in isoform_data.iterrows():
        # Strip suffixes from isoform IDs before comparison
        reviewed_isoforms = strip_suffix(entry.get("reviewed_isoforms", ""))
        unreviewed_isoforms = strip_suffix(entry.get("unreviewed_isoforms", ""))
        
        # Check if the UniProt ID matches any isoform (stripped version)
        if uniprot == reviewed_isoforms or uniprot == unreviewed_isoforms:
            # If a match is found, add the uniprotkb_canonical_ac to the result
            uniprotkb_ac = strip_suffix(entry.get("uniprotkb_canonical_ac"))
            
            # Check if the UniProt ID matches the canonical version
            if uniprot == uniprotkb_ac:
                is_canonical = "yes"
            
            # Store the result with canonical status
            result[ensp] = {"uniprotkb_canonical_ac": uniprotkb_ac, "canonical": is_canonical}
            break  # Exit inner loop once the first match is found

# Write the result to a JSON file
with open(canonical_toy_output_path, "w") as json_file:
    json.dump(result, json_file, indent=4)
