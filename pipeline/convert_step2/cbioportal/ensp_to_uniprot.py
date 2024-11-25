import requests
import time
import pandas as pd

"""
# Step 1: Load ENSEMBL protein IDs from a CSV file
file_path = "/data/shared/biomuta/generated/datasets/2024_10_22/mapping_ids/chr_pos_to_ensp_toy.csv"
df = pd.read_csv(file_path)
ensp_ids = df['ENSP'].tolist()  # 'ENSP' is the column name for ENSP IDs
ids_str = ",".join(ensp_ids)
"""

# Step 2: Submit the ID mapping job
url = "https://rest.uniprot.org/idmapping/run"
params = {
    "from": "Ensembl_Protein",
    "to": "UniProtKB_AC-ID",
    "ids": "ENSP00000359446"#ids_str
}

response = requests.post(url, data=params)
if response.ok:
    job_id = response.json()["jobId"]
    print(f"Job ID: {job_id}")

    # Step 3: Poll for job completion
    result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    while True:
        result_response = requests.get(result_url)
        if result_response.ok:
            result_data = result_response.json()
            if result_data.get("jobStatus") == "FINISHED":
                print("Mapping job completed.")
                break
        else:
            print("Waiting for job completion...")
        time.sleep(2)  # Poll every 2 seconds

    # Step 4: Retrieve and process the results
    mapped_results = result_data.get("results", [])
    ensp_to_uniprot = {entry['from']: entry['to']['primaryAccession'] for entry in mapped_results}
    canonical_accessions = [entry['to']['primaryAccession'] for entry in mapped_results]
    print("Canonical UniProt Accessions:", canonical_accessions)
    """
    # Step 5: Map UniProt accessions back to the original DataFrame
    df['UniProt_Accession'] = df['ENSP_ID'].map(ensp_to_uniprot)  # Add the new column

    # Step 6: Save the updated DataFrame to a new CSV file
    updated_file_path = "/data/shared/biomuta/generated/datasets/2024_10_22/mapping_ids/ensp_to_uniprot_toy.csv"  # Updated file
    df.to_csv(updated_file_path, index=False)
    print(f"Updated CSV file saved to {updated_file_path}")
    """
else:
    print("Failed to submit the job.")
