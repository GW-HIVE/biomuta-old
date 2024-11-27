import json
import logging
import os
import sys
from pathlib import Path
from typing import Optional

logging.basicConfig(filename="cancer_mapping.log",
                    filemode='a',
                    format='%(asctime)s %(levelname)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO)

# Logging levels
# 1. error
# 2. warning
# 3. info
# 4. debug

# Add this to .gitignore
# *.log

logging.info("Logger started ----------------------")
# Put this code into a utils file, define frequently used variables in a util python script

#sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
#from utils import ROOT_DIR #__init__.py, ROOT_DIR is a global var

# Define paths

# Get the directory of this script
script_dir = Path(__file__).resolve().parent
# Navigate to config.json location relative to script
config_dir = script_dir.parent.parent
# Load config
with open(config_dir/'config.json') as config_file:
    config = json.load(config_file)
# Access paths from config
mapping_dir = Path(config["relevant_paths"]["mapping"])
doid_mapping = mapping_dir / "combined_do_mapping.json"
fallback_do_map = mapping_dir / "fallback_cbio_doid_mapping.json"

# Input and output file names
# Get the latest directory
directory_path = Path(config["relevant_paths"]["generated_datasets"])
latest_dir = max([d for d in os.listdir(directory_path) if os.path.isdir(os.path.join(directory_path, d))], key=lambda d: os.path.getctime(os.path.join(directory_path, d)))
latest_dir = Path(directory_path) / latest_dir
def ask_confirmation(prompt):
    while True:
        user_input = input(f"{prompt} (y/n): ").strip().lower()
        if user_input == 'y':
            return True
        elif user_input == 'n':
            return False
        else:
            print(f"Invalid input. Please enter 'y' for yes or 'n' for no.")
if ask_confirmation(f"The latest created directory is: {latest_dir}. Proceed?"):
    input_file = Path(latest_dir) / "unique_cancer_names.txt"
    cancer_types_with_do = Path(latest_dir) / "cancer_types_with_do.json"
    cancer_type_per_study = Path(latest_dir) / "cancer_type_per_study.json"
    study_ids_with_do = Path(latest_dir) / "study_ids_with_do.json"
    print(f"Using {latest_dir}/unique_cancer_names.txt and writing out to {latest_dir}/cancer_types_with_do.json")
else:
    sys.exit("Aborted by user.")

# Open and load the primary and fallback mapping files
with open(doid_mapping, "r") as f:
    map1 = json.load(f)
with open(fallback_do_map, "r") as map2_f:
    map2 = json.load(map2_f)

def check_do_mapping(cancer_type: str) -> str:
    cancer_type = cancer_type.lower()

    # Check in primary (civic and cosmic) mapping file first
    for source in ["civic", "cosmic"]:
        for entry in map1[source]:
            if entry["cancer_name"] == cancer_type:
                logging.info(f"{cancer_type} matched")
                return entry["do_term"]
    
    # If no match in primary, check fallback mapping
    matches: list[str] = []
    keywords = map2.keys()
    for keyword in keywords:
        if keyword == cancer_type:
            matches = [map2[keyword]]
            break
        elif keyword in cancer_type.split():
            matches.append(map2[keyword])
    if len(matches) > 1:
        logging.warning(f"Cancer type: {cancer_type} matched on {matches}")
        return "Too many cancer type matches."
    elif not matches:
        logging.warning(f"No matches for {cancer_type}")
        return "NA"
    else:
        return matches[0]

# Read the input file and generate the JSON objects
cancer_list = []
with open(input_file, "r") as f:
    for line in f:
        cancer_type = line.strip()
        if cancer_type:  # Skip empty lines
            do_name = check_do_mapping(cancer_type)
            cancer_list.append({"cancerType": cancer_type, "do_name": do_name})

# Write the JSON data to the output file
with open(cancer_types_with_do, "w") as f:
    json.dump(cancer_list, f, indent=2)

logging.info(f"JSON file generated: {cancer_types_with_do}")

"""
if keyword in cancer_type:
  match = map[keyword]
  returning match

  mapped_cancer_term
  full_map[resource][cancer_type] = mapped_cancer_term
  json.dump(open("path", "w"), full_map)
"""
#make a Venn diagram for EGFR with mutations and see which mutations in the gene are in BioMuta and which are in cBio and how many are overlapping
#how many new slim did I find?

# Load cancer_types_with_do
with open(cancer_types_with_do, 'r') as f:
    cancer_do_data = json.load(f)

# Load cancer_type_per_study
with open(cancer_type_per_study, 'r') as f:
    cancer_study_data = json.load(f)

# Create a dictionary to map 'cancerType' to 'do_name' from cancer_types_with_do
do_name_map = {entry['cancerType']: entry['do_name'] for entry in cancer_do_data}

# Create the result list with "studyId" and corresponding "do_name"
result = [
    {
        "studyId": entry['studyId'],
        "do_name": do_name_map.get(entry['cancerType'])
    }
    for entry in cancer_study_data if entry['cancerType'] in do_name_map
]

# Save the result to a new JSON file
with open(study_ids_with_do, 'w') as f:
    json.dump(result, f, indent=2)

logging.info(f"Output JSON file {study_ids_with_do} created with 'studyId' and 'do_name'")
