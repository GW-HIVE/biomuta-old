#!/bin/bash

# Input and output file paths
input_csv="/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/chr_pos_to_ensp.csv"      # Input CSV file
output_json="/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/ensp_to_uniprot_mappings.json"   # Output JSON file

batch_size=5000            # Number of ENSP IDs per batch (adjustable)

# Create a temporary file for unique ENSP IDs
temp_file=$(mktemp)

# Extract unique ENSP IDs from the input CSV (excluding the header)
awk -F, 'NR > 1 {print $6}' "$input_csv" | sort -u > "$temp_file"

# Create an associative array to store mappings
declare -A mappings

# Ensure cleanup of temporary files on exit
trap 'rm -f "$temp_file" "$batch_file" batch_*' EXIT

# Process ENSP IDs in batches
batch_file=$(mktemp)  # Temporary file for batch data
split -l "$batch_size" "$temp_file" batch_

for batch_chunk in batch_*; do
    # Create a comma-separated list of ENSP IDs for the batch
    paste -sd, "$batch_chunk" > "$batch_file"

    # Submit the batch request to the UniProt API
    response=$(curl --silent --request POST 'https://rest.uniprot.org/idmapping/run' \
        --form "ids=$(cat $batch_file)" \
        --form 'from="Ensembl_Protein"' \
        --form 'to="UniProtKB"')
    jobId=$(echo "$response" | jq -r '.jobId')

    if [ -z "$jobId" ]; then
        echo "Failed to retrieve jobId for batch: $(cat $batch_chunk)"
        continue
    fi

    echo "Processing batch: $(basename "$batch_chunk")"

    # Poll the job status
    while true; do
        status_response=$(curl --silent "https://rest.uniprot.org/idmapping/status/$jobId")
        status=$(echo "$status_response" | jq -r '.jobStatus')
        
        if [ "$status" == "FINISHED" ]; then
            break
        elif [ "$status" == "FAILED" ]; then
            echo "Job failed for batch: $(cat $batch_chunk)"
            continue 2
        else
            sleep 5
        fi
    done

    # Fetch and parse results
    result=$(curl --silent "https://rest.uniprot.org/idmapping/uniprotkb/results/$jobId")
    if [ -z "$result" ]; then
        echo "No results for batch: $(cat $batch_chunk)"
        continue
    fi

    # Process API results and populate mappings
    echo "$result" | jq -c '.results[]' > temp_results.json

    while read -r record; do
        ensp_id=$(echo "$record" | jq -r '.from')
        primaryAccession=$(echo "$record" | jq -r '.to.primaryAccession')
        mappings["$ensp_id"]=$primaryAccession
    done < temp_results.json

    rm -f temp_results.json  # Clean up
done

# Write the mappings to a JSON file in the desired format
echo "{" > "$output_json"
for key in "${!mappings[@]}"; do
    echo "  \"$key\": \"${mappings[$key]}\"," >> "$output_json"
done
# Remove the last comma from the last entry
sed -i '$ s/,$//' "$output_json"
echo "}" >> "$output_json"

echo "Processing complete. Output saved to $output_json"
