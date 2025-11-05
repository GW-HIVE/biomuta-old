import csv

bm_path = '/data/shared/repos/biomuta-old/generated_datasets/compiled/O14497-1.csv'
clinical_path = '/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/clinical_information.csv'
output_path = '/data/shared/repos/biomuta-old/generated_datasets/compiled/O14497-1_clinical.csv'

# Read first column from each csv
def read_first_column(file_path):
    with open(file_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        return [row[0] for row in reader if row]
    
bm_fields = set(read_first_column(bm_path))
print(f"bm_fields length: {len(bm_fields)}")
clinical_fields = set(read_first_column(clinical_path))
print(f"clinical_fields length: {len(clinical_fields)}")

# Find matches
common_fields = bm_fields.intersection(clinical_fields)
print(f"common_fields length: {len(common_fields)}")
"""
with open(output_path, 'w') as f:
    for field in common_fields:
        f.write(f"{field}\n")
"""