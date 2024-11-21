import os
import json
import glob
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
from utils import ROOT_DIR
from utils.config import get_config

def process_json_to_bed(input_directory, output_bed_file):
    buffer = []
    buffer_size = 1000  # Adjust based on memory availability
    file_count = 0
    error_count = 0

    with open(output_bed_file, 'w') as bed_file:
        for json_file_path in glob.glob(os.path.join(input_directory, '*.json')):
            file_count += 1
            try:
                with open(json_file_path, 'r') as json_file:
                    # Load JSON data
                    data = json.load(json_file)
                    
                    # Process each JSON object in the file
                    for record in data:
                        try:
                            # Extract chromosome, start position, and end position
                            chr_ = record['chr']
                            # Replace chr23 with chrX
                            if chr_ == '23':
                                chr_ = 'X'
                            # Ensure chromosome name starts with 'chr' (if not already)
                            if not chr_.startswith('chr'):
                                chr_ = 'chr' + chr_
                            start_pos = record['startPosition'] - 1  # Convert to 0-based for BED
                            end_pos = record['endPosition']

                            # Append line to buffer
                            buffer.append(f"{chr_}\t{start_pos}\t{end_pos}\n")

                            # Write buffer to file when it reaches the specified size
                            if len(buffer) >= buffer_size:
                                bed_file.writelines(buffer)
                                buffer.clear()

                        except KeyError as e:
                            error_count += 1
                            print(f"Error: Missing key {e} in record from file {json_file_path}")
                        except Exception as e:
                            error_count += 1
                            print(f"Error: An error occurred while processing record in file {json_file_path}: {e}")

            except json.JSONDecodeError:
                error_count += 1
                print(f"Error: Could not decode JSON from file {json_file_path}")
            except Exception as e:
                error_count += 1
                print(f"Error: An error occurred while reading file {json_file_path}: {e}")

            # Output progress every 10 files processed
            if file_count % 10 == 0:
                print(f"Processed {file_count} files...")

        # Write any remaining data in the buffer
        if buffer:
            bed_file.writelines(buffer)

    print(f"Processing complete. Total files processed: {file_count}, errors encountered: {error_count}")

# Run the function
config_obj = get_config()
dl_dir = Path(config_obj["relevant_paths"]["downloads"])
out_dir = Path(config_obj["relevant_paths"]["generated_datasets"])
input_directory = dl_dir / 'cbioportal' / '2024_10_21' / 'mutations'  # Write a util to get latest dir
output_bed_file = out_dir / '2024_10_22' / 'hg19positions.bed' #Write a util to get latest dir
process_json_to_bed(input_directory, output_bed_file)
