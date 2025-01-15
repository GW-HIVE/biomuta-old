import pandas as pd

# Load the biomuta.csv file
biomuta_df = pd.read_csv("/data/shared/biomuta/downloads/biomuta.csv")

# Load the human_protein_masterlist.csv file
human_protein_df = pd.read_csv("/data/shared/repos/biomuta-old/downloads/glygen/human_protein_masterlist.csv")

# Get unique uniprotkb_canonical_ac values from human_protein_masterlist.csv
masterlist_uniprot_set = set(human_protein_df["uniprotkb_canonical_ac"].dropna())
reviewed_isoforms_set = set(human_protein_df["reviewed_isoforms"].dropna())
unreviewed_isoforms_set = set(human_protein_df["unreviewed_isoforms"].dropna())

# Convert biomuta uniprotkb_canonical_ac column to a set to remove duplicates
unique_biomuta_uniprot_set = set(biomuta_df['uniprotkb_canonical_ac'].dropna())
print(f"Unique UniProt ACs in biomuta.csv: {len(unique_biomuta_uniprot_set)}")

# Check if each uniprotkb_canonical_ac in biomuta.csv is in the masterlist
def check_isoforms(uniprotkb_ac):
    if uniprotkb_ac in masterlist_uniprot_set:
        return 'Found in masterlist'
    elif uniprotkb_ac in reviewed_isoforms_set:
        return 'Found in reviewed_isoforms'
    elif uniprotkb_ac in unreviewed_isoforms_set:
        return 'Found in unreviewed_isoforms'
    else:
        return 'Not found in either'

# Apply the check_isoforms function to each unique uniprotkb_canonical_ac
biomuta_results = {
    'uniprotkb_canonical_ac': list(unique_biomuta_uniprot_set),
    'isoform_status': [check_isoforms(ac) for ac in unique_biomuta_uniprot_set]
}

# Convert the results into a DataFrame
biomuta_results_df = pd.DataFrame(biomuta_results)

# Print the counts of found and not found
found_count = biomuta_results_df['isoform_status'].apply(lambda x: 'Found in masterlist' in x).sum()
not_found_count = len(biomuta_results_df) - found_count

print(f"Found in masterlist: {found_count}")
print(f"Not found in masterlist: {not_found_count}")

# Save the result to a new CSV file
biomuta_results_df.to_csv("/data/shared/biomuta/downloads/biomuta_with_masterlist_check_and_isoforms.csv", index=False)
