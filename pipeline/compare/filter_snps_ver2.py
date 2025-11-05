import json
import glob

# Directory containing mutation JSON files
input_dir = '/data/shared/repos/biomuta-old/downloads/cbioportal/2024_10_21/mutations'
output_file = '/data/shared/repos/biomuta-old/downloads/cbioportal/2024_10_21/filtered_mutations/filtered_snps_ver2.json'

# Write an empty list to initialize the JSON output file
with open(output_file, 'w') as outfile:
    outfile.write('[')

# Collect all files in the directory
json_files = glob.glob(f"{input_dir}/*.json")
first_entry = True

for file_path in json_files:
    with open(file_path, 'r') as file:
        data = json.load(file)
        
        for mutation in data:
            if mutation.get('startPosition') == mutation.get('endPosition'):
                snp_entry = {
                    "uniquePatientKey": mutation.get('uniquePatientKey'),
                    "patientId": mutation.get('patientId'),
                    "entrezGeneId": mutation.get('entrezGeneId'),
                    "studyId": mutation.get('studyId'),
                    "chr": mutation.get('chr'),
                    "startPosition": mutation.get('startPosition'),
                    "endPosition": mutation.get('endPosition'),
                    "proteinChange": mutation.get('proteinChange'),
                    "refseqMrnaId": mutation.get('refseqMrnaId')
                }
                
                # Append comma before each entry except the first one
                if not first_entry:
                    with open(output_file, 'a') as outfile:
                        outfile.write(',\n')
                else:
                    first_entry = False
                
                # Append the entry to the output file
                with open(output_file, 'a') as outfile:
                    json.dump(snp_entry, outfile)

# Close the JSON array properly
with open(output_file, 'a') as outfile:
    outfile.write(']')
