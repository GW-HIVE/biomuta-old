import pandas as pd

# Load the CSV files
biomuta = pd.read_csv("/data/shared/biomuta/generated/stats/full_tables/biomuta_canonical.csv")
venn = pd.read_csv("/data/shared/biomuta/generated/stats/full_tables/venn_4_components.tsv", sep="\t")

# Filter rows in venn_4_components.tsv with "only_in_biomuta"
venn_only_biomuta = venn[venn['Source'] == 'only_in_biomuta']
print(len(venn_only_biomuta))

# Remove duplicates from biomuta based on the relevant columns
biomuta_unique = biomuta.drop_duplicates(subset=['uniprotkb_canonical_ac', 'aa_pos', 'ref_aa', 'alt_aa'])

# Perform an inner merge between biomuta and venn data
merged = pd.merge(
    biomuta_unique,
    venn_only_biomuta,
    left_on=['uniprotkb_canonical_ac', 'aa_pos', 'ref_aa', 'alt_aa'],
    right_on=['uniprotkb_canonical_ac', 'aa_pos', 'ref_aa', 'alt_aa']
)

# Count occurrences of unique values in the "source" column
source_counts = merged['source'].value_counts()

# Display the results
print(source_counts)