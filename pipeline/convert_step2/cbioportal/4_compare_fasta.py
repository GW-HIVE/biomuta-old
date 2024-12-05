import requests
import json
from pathlib import Path

# Load config.json
config_path = Path(__file__).resolve().parent.parent.parent / "config.json"
with open(config_path, "r") as config_file:
    config = json.load(config_file)

# Retrieve  paths from config.json
repos_base = Path(config["relevant_paths"]["repos_generated_datasets"])
ensembl_uniprot_map_path = repos_base / "2024_10_22/mapping_ids/canonical_toy.json"

# Load your JSON file containing ENSEMBL to UniProt mappings
with open(ensembl_uniprot_map_path, "r") as file:
    ensembl_uniprot_map = json.load(file)

def fetch_ensembl_sequence(ensembl_id):  # Fetch ENSEMBL sequence
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    else:
        print(f"Failed to fetch ENSEMBL sequence for {ensembl_id}")
        return None

def fetch_uniprot_sequence(uniprot_id):  # Fetch UniProt sequence
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.split('\n', 1)[1].replace('\n', '')
    else:
        print(f"Failed to fetch UniProt sequence for {uniprot_id}")
        return None

# Compare sequences
for ensembl_id, uniprot_id in ensembl_uniprot_map.items():
    ensembl_sequence = fetch_ensembl_sequence(ensembl_id)
    uniprot_sequence = fetch_uniprot_sequence(uniprot_id)
    
    if ensembl_sequence and uniprot_sequence:
        if ensembl_sequence == uniprot_sequence:
            print(f"Sequences match for {ensembl_id} and {uniprot_id}")
        else:
            print(f"Sequences do not match for {ensembl_id} and {uniprot_id}")
    else:
        print(f"Could not compare sequences for {ensembl_id} and {uniprot_id}")
