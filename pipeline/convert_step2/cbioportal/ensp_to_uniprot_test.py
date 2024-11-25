import requests

def map_ensembl_to_uniprot(ensembl_ids):
    # Set up your API endpoint and headers if needed
    url = "https://rest.uniprot.org/idmapping"
    headers = {"Content-Type": "application/json"}

    # Send request to the mapping API
    response = requests.post(url, headers=headers, json={"from": "Ensembl_Protein", "to": "UniProtKB", "ids": ensembl_ids})
    data = response.json()

    # Extract only the primary accession IDs from results
    primary_accessions = [result['to']['primaryAccession'] for result in data['results'] if 'to' in result and 'primaryAccession' in result['to']]
    
    return primary_accessions

# Example usage:
ensembl_ids = ["ENSP00000359446"]  # List of ENSEMBL IDs
primary_accessions = map_ensembl_to_uniprot(ensembl_ids)
print(primary_accessions)
