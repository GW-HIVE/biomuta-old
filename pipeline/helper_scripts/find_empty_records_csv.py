import csv

# Specify the CSV file and the column name you want to check
csv_file = "/data/shared/biomuta/generated/datasets/2024_10_22/mapping_ids/chr_pos_to_ensp_old.csv"
column_name = "prot_change"

# Open the CSV file and read it
with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        # Check if the value in the specified column is empty
        if not row[column_name].strip():  # .strip() removes leading/trailing whitespaces
            print(f"Empty record found in row: {row}")

