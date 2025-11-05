import pandas as pd

# Read the data from the file
file_path = "/data/shared/biomuta/generated/stats/full_tables/uniprot_ac_cancer_type.tsv"
df = pd.read_csv(file_path, sep="\t")

# Helper function to extract DOIDs
def extract_doids(cancer_type):
    return set(part.split(" /")[0] for part in cancer_type.split(" & "))

# Add a column with extracted DOID sets
df["DOID_Set"] = df["do_name"].apply(extract_doids)

# Group by the first column and find partial matches
results = []
for protein, group in df.groupby("uniprotkb_canonical_ac"):
    doid_sets = group["DOID_Set"].tolist()
    if len(doid_sets) > 1:  # Check only if there are multiple rows for a protein
        for i in range(len(doid_sets)):
            for j in range(i + 1, len(doid_sets)):
                if doid_sets[i] & doid_sets[j]:  # Check intersection of DOID sets
                    results.append(group.iloc[[i, j]])

# Combine results
if results:
    matched_rows = pd.concat(results)
    matched_rows.to_csv("/data/shared/biomuta/generated/stats/full_tables/partial_matches.tsv", sep="\t", index=False)
    print("Partial matches saved to partial_matches.tsv")
else:
    print("No partial matches found.")

