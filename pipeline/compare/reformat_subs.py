import csv
import re

# Input and output file paths
input_file = '/data/shared/biomuta/generated/datasets/2024_10_22/cbio_egfr_sub.txt'
output_file = '/data/shared/biomuta/generated/datasets/2024_10_22/cbio_egfr_sub_reshaped.csv'

# Regular expression pattern to match lines in your specified format
pattern = r"([A-Z])(\d+)([A-Z\*])"

# List to hold the parsed data for sorting
data = []

with open(input_file, 'r') as infile:
    for line in infile:
        match = re.match(pattern, line.strip())
        if match:
            # Extracting groups from regex
            first_letter = match.group(1)
            number = match.group(2)
            second_letter = match.group(3)
            # Append to data list
            data.append([int(number), first_letter, second_letter])

# Sort data by the first column (numeric part)
data.sort(key=lambda x: x[0])

# Write the sorted data to the CSV file
with open(output_file, 'w', newline='') as outfile:
    csv_writer = csv.writer(outfile)
    csv_writer.writerows(data)