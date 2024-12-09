import csv
import json
import os
import pandas as pd
from collections import Counter

# Function to count and print the occurrences of each record
def count_records(json_files, record_name):
    record_counter = Counter()  # To store the count of each unique record

    for json_file in json_files:
        with open(json_file, 'r') as f:
            data = json.load(f)  # Load the JSON data

            # Iterate over each record and count record
            for record in data:
                if record_name in record:
                    record_counter[record[record_name]] += 1  # Increment count for the record

    # Print the count of each record
    print(f"{record_name} counts:")
    for record, count in record_counter.items():
        print(f"{record}: {count}")

'''
# Example usage
# Get a list of all JSON files in the directory
json_files = [os.path.join(json_directory, f) for f in os.listdir(json_directory) if f.endswith('.json')]
# Call the function to count variant types
count_records(json_files, 'chr') #Change to pass record_name as a command line argument
'''

def count_unique_values(file_path, column_name, output_file):
    """
    Count occurrences of each unique value in a specified column of a CSV file
    and write the results to a TSV file.

    Parameters:
    - file_path (str): Path to the CSV file.
    - column_name (str): The name of the column to analyze.
    - output_file (str): Path to the output TSV file.
    """
    # Load the CSV file into a DataFrame
    df = pd.read_csv(file_path)

    # Count the occurrences of each unique value in the specified column
    value_counts = df[column_name].value_counts()

    # Write the results to a TSV file
    value_counts.to_csv(output_file, sep='\t', header=['Count'], index_label=column_name)
'''
# Example usage:
count_unique_values('/data/shared/biomuta/downloads/biomuta.csv', 'do_name', '/data/shared/biomuta/generated/stats/sites_per_do_bm.tsv')
count_unique_values('/data/shared/biomuta/generated/datasets/2024_10_22/mapping_ids/chr_pos_to_ensp_old_do_name.csv', 'do_name', '/data/shared/biomuta/generated/stats/sites_per_do_cbio.tsv')
'''


def compare_tsv_and_output(file1, file2, output_file):
    # Load the first file into a dictionary
    data1 = {}
    with open(file1, 'r') as f1:
        reader = csv.reader(f1, delimiter='\t')
        next(reader)  # Skip the header line
        for row in reader:
            data1[row[0]] = int(row[1])

    # Load the second file into a dictionary
    data2 = {}
    with open(file2, 'r') as f2:
        reader = csv.reader(f2, delimiter='\t')
        next(reader)  # Skip the header line
        for row in reader:
            data2[row[0]] = int(row[1])

    # Open the output file to write the comparison results
    with open(output_file, 'w', newline='') as f_out:
        writer = csv.writer(f_out, delimiter='\t')
        writer.writerow(['do_name', 'cbio_site_count', 'biomuta_site_count'])

        # Compare the data and write the results
        all_keys = set(data1.keys()).union(data2.keys())
        for key in all_keys:
            value1 = data1.get(key, 'N/A')
            value2 = data2.get(key, 'N/A')
            writer.writerow([key, value1, value2])

# Example usage:
compare_tsv_and_output('/data/shared/biomuta/generated/stats/sites_per_do_cbio.tsv', '/data/shared/biomuta/generated/stats/sites_per_do_bm.tsv', '/data/shared/biomuta/generated/stats/sites_per_do_bm_cbio.tsv')


def find_matching_json_records(directory, chr_value, entrezGeneId_value):
    """
    Find JSON records in json files in the specified directory that have 'chr' as chr_value
    and 'entrezGeneId' as entrezGeneId_value.

    Args:
        directory (str): The path to the directory containing the JSON files.
        chr_value (str): The value of the 'chr' field to match, should be a string.
        entrezGeneId_value (int): The value of the 'entrezGeneId' field to match, should be an integer.

    Returns:
        list: A list of matching json records.
    """
    # Ensure chr_value is a string and entrezGeneId_value is an integer
    chr_value = str(chr_value)
    entrezGeneId_value = int(entrezGeneId_value)

    matching_records = []

    # Iterate over each file in the directory
    for filename in os.listdir(directory):
        if filename.endswith('.json'):
            file_path = os.path.join(directory, filename)
            
            # Open and read the JSON file
            with open(file_path, 'r') as f:
                try:
                    data = json.load(f)

                    # Iterate over each dictionary in the list
                    for item in data:
                        if isinstance(item, dict):
                            if item.get('chr') == chr_value and item.get('entrezGeneId') == entrezGeneId_value:
                                matching_records.append(item)

                except json.JSONDecodeError:
                    print(f"Skipping invalid JSON file: {file_path}")

    return matching_records
'''
# Example usage
json_directory = '/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations'
chr_value = 15  # This will be converted to a string
entrezGeneId_value = 28472  # This should be an integer
matching_records = find_matching_json_records(json_directory, chr_value, entrezGeneId_value)
print(f"Matching entries: {matching_records}")
'''

def find_records(filename, chr_value, entrezGeneId_value):
    # Ensure chr_value is a string and entrezGeneId_value is an integer
    chr_value = str(chr_value)
    entrezGeneId_value = int(entrezGeneId_value)

    records = []

    with open(filename, 'r') as f:
        data = json.load(f)

        for item in data:
            if isinstance(item, dict):
                if item.get('chr') == chr_value and item.get('entrezGeneId') == entrezGeneId_value:
                    records.append(item)

    return records