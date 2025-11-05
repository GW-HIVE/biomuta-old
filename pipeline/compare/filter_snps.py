import json
import glob

# Directory containing mutation JSON files
input_dir = '/data/shared/repos/biomuta-old/downloads/cbioportal/2024_10_21/mutations'
output_file = '/data/shared/repos/biomuta-old/downloads/cbioportal/2024_10_21/filtered_mutations/filtered_snps.json'

# Collect all files in the input directory
json_files = glob.glob(f"{input_dir}/*.json")

# Initialize an empty list to store SNPs
snps = []

# Process each file
for file_path in json_files:
    with open(file_path, 'r') as file:
        data = json.load(file)

        # Iterate over each mutation in the file
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
                    "ncbiBuild": mutation.get('ncbiBuild'),
                    "proteinChange": mutation.get('proteinChange'),
                    "refseqMrnaId": mutation.get('refseqMrnaId')
                }
                snps.append(snp_entry)

# Write the filtered SNP data to a JSON output file
with open(output_file, 'w') as outfile:
    json.dump(snps, outfile, indent=4)