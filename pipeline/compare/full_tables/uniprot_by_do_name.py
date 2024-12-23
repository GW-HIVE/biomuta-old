import pandas as pd

# Load the CSV files
cbio = pd.read_csv("/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/final_table.csv")
biomuta = pd.read_csv("/data/shared/repos/biomuta-old/downloads/biomuta.csv")

# Perform the comparisons
cbio_set = set(cbio["uniprotkb_canonical_ac"])
biomuta_set = set(biomuta["uniprotkb_canonical_ac"])

in_both = cbio_set & biomuta_set
only_in_cbio = cbio_set - biomuta_set
only_in_biomuta = biomuta_set - cbio_set

# Create a DataFrame for the results
results = []

# Add "in_both" entries
for uniprot in in_both:
    results.append({"uniprotkb_canonical_ac": uniprot, "category": "in_both"})

# Add "only_in_cbio" entries
for uniprot in only_in_cbio:
    results.append({"uniprotkb_canonical_ac": uniprot, "category": "only_in_cbio"})

# Add "only_in_biomuta" entries
for uniprot in only_in_biomuta:
    results.append({"uniprotkb_canonical_ac": uniprot, "category": "only_in_biomuta"})

# Convert to a DataFrame
results_df = pd.DataFrame(results)

# Save to a TSV file
results_df.to_csv("uniprot_ac_by_source.tsv", sep="\t", index=False)

print("Comparison results saved to 'uniprot_ac_by_source.tsv'")
