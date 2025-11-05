import requests
import json

# Load your JSON file containing ENSEMBL to UniProt mappings
with open("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/canonical.json") as file:
    ensembl_uniprot_map = json.load(file)

def fetch_ensembl_sequence(ensembl_id): #it works
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/plain"
    #https://rest.ensembl.org/sequence/id/{ensembl_id}?content-type=text/x-fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    else:
        print(f"Failed to fetch ENSEMBL sequence for {ensembl_id}")
        return None

def fetch_uniprot_sequence(uniprot_id): #it works
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.split('\n', 1)[1].replace('\n', '')
    else:
        print(f"Failed to fetch UniProt sequence for {uniprot_id}")
        return None

# Compare sequences and write the results
with open("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/fasta_results.txt", "w") as output_file:
    for ensembl_id, uniprot_id in ensembl_uniprot_map.items():
        ensembl_sequence = fetch_ensembl_sequence(ensembl_id)
        uniprot_sequence = fetch_uniprot_sequence(uniprot_id)

        if ensembl_sequence and uniprot_sequence:
            if ensembl_sequence == uniprot_sequence:
                output_file.write(f"Sequences match for {ensembl_id} and {uniprot_id}\n")
            else:
                output_file.write(f"Sequences do not match for {ensembl_id} and {uniprot_id}\n")
        else:
            output_file.write(f"Could not compare sequences for {ensembl_id} and {uniprot_id}\n")
