import requests

def fetch_uniprot_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.split('\n', 1)[1].replace('\n', '')
    else:
        print(f"Failed to fetch UniProt sequence for {uniprot_id}")
        return None
    
print(fetch_uniprot_sequence('Q9NXB0'))