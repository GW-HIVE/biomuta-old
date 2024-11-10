import csv
import json

# Define dictionaries for each file's data
data = {
    "civic": [],
    "cosmic": []
}

# Track seen disease names for each file to remove duplicates
seen_names = {
    "civic": set(),
    "cosmic": set()
}

# Process civic and cosmic files
civic_cosmic = [
    ("/data/shared/biomuta/pipeline/convert_step2/mapping/civic_doid_mapping.csv", "civic"),
    ("/data/shared/biomuta/pipeline/convert_step2/mapping/cosmic_doid_mapping.csv", "cosmic")
]

for input_file, source_key in civic_cosmic:
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=',')
        
        # Skip the header
        next(reader)
        
        for row in reader:
            # Process the first column: make lowercase, replace underscores with spaces, and remove the part in parentheses
            first_column = row[0].split('(')[0].replace('_', ' ').strip().lower()
            
            # Only add to data if this cancer_name has not been seen before in this source
            if first_column not in seen_names[source_key]:
                data[source_key].append({
                    "cancer_name": first_column,
                    "do_term": row[1]
                })
                # Mark this cancer_name as seen
                seen_names[source_key].add(first_column)

# Save the data dictionary to a JSON file
with open("/data/shared/biomuta/pipeline/convert_step2/mapping/combined_do_mapping.json", 'w') as outfile:
    json.dump(data, outfile, indent=4)
