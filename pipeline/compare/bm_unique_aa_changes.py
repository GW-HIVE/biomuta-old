import csv

# Open the CSV file
input_file = '/data/shared/biomuta/downloads/egfr.csv'  # replace with your CSV file path
output_file = '/data/shared/biomuta/generated/datasets/2024_10_22/bm_egfr_unique_combinations.txt'  # output file to write the results

# Initialize a set to store unique combinations
unique_combinations = set()

# Open the CSV and process it
with open(input_file, mode='r') as file:
    reader = csv.reader(file)
    for row in reader:
        # Extract values from columns 7, 8, and 9 (indexes 6, 7, and 8)
        combination = (row[6], row[7], row[8])
        unique_combinations.add(combination)

# Write the unique combinations to a text file
with open(output_file, mode='w') as file:
    for combination in unique_combinations:
        file.write(f"{combination[0]}, {combination[1]}, {combination[2]}\n")

print(f"Unique combinations have been written to {output_file}")
