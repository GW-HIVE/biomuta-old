import pandas as pd

# Load the CSV files
biomuta = pd.read_csv("/data/shared/biomuta/downloads/biomuta.csv", quotechar='"')
masterlist_check = pd.read_csv("/data/shared/biomuta/downloads/biomuta_with_masterlist_check_and_isoforms.csv")

# Filter for "Found in masterlist"
found_in_masterlist = masterlist_check[masterlist_check["isoform_status"] == "Found in masterlist"]

# Extract uniprotkb_canonical_ac values
masterlist_set = set(found_in_masterlist["uniprotkb_canonical_ac"])

# Filter biomuta.csv rows
biomuta_filtered = biomuta[biomuta["uniprotkb_canonical_ac"].isin(masterlist_set)]

# Save filtered data to a new file
biomuta_filtered.to_csv("/data/shared/biomuta/generated/stats/full_tables/biomuta_canonical.csv", index=False)
