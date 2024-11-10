import ijson
import csv
import traceback
import os
import glob

# Specify the directory containing the JSON files
json_directory = "/data/shared/pipelines/cbioportal/mutations"

# Get all JSON files in the specified directory
json_files = glob.glob(os.path.join(json_directory, "*.json"))

# List to keep track of files that encounter errors
error_files = []

# Outer loop: Process each JSON file
for json_file in json_files:
    try:
        print(f"Starting file {json_file}")

        # Create the output TSV file
        output_filename = f"{os.path.splitext(json_file)[0]}_data.tsv"
        with open(output_filename, "w") as tsvfile:
            
            # Open the JSON file for reading
            with open(json_file, "r") as file:
                
                # Initialize the TSV writer
                tsv_writer = None
                
                # Loop through each record in the JSON file using ijson
                for idx, record in enumerate(ijson.items(file, "item")):
                    
                    # For the first record, extract the headers and write them
                    if idx == 0:
                        headers = set(record.keys())  # Extract headers from the first record
                        tsv_writer = csv.DictWriter(tsvfile, delimiter="\t", fieldnames=headers)
                        tsv_writer.writeheader()  # Write the header to the TSV file
                    
                    # Write the record to the TSV file
                    try:
                        tsv_writer.writerow(record)
                    except Exception as e:
                        # Handle any errors during writing
                        traceback.print_exc()
                        print(f"Failed to write row {idx} from file {json_file}")
                        error_files.append(json_file)
                        break

    except Exception as e:
        # Handle any errors in opening/processing the file
        traceback.print_exc()
        print(f"Error processing file {json_file}")
        error_files.append(json_file)

# After processing all files, print the list of files that encountered errors
print(f"Files with errors: {error_files}")