import json
import os
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

# Directory containing your JSON files (adjust as needed)
json_directory = '/data/shared/biomuta/downloads/cbioportal/2024_10_21/mutations'

# Get a list of all JSON files in the directory
json_files = [os.path.join(json_directory, f) for f in os.listdir(json_directory) if f.endswith('.json')]

# Call the function to count variant types
count_records(json_files, 'chr') #Change to pass record_name as a command line argument
