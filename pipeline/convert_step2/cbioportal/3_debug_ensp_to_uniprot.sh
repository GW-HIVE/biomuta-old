#!/bin/bash

# Successfully created batches of ESNP IDs of size 5000
# Testing API calls
# API only processes 25 IDs at a time for some reason
# Reducing batch size to 25
# It works for the first batch, testing all batches
# An error appears at batch 879
# Stopped there
# Mapping is going to be done using human_protein_transcriptlocus.csv from data.glygen.org/GLY_000135 (see 3_ensp_to_uniprot.py) to avoid using the API.
# 87,679 ENSP IDs were mapped, 23,517 ENSP IDs remain unmapped (see if unmapped IDs are non-canonical = process separately).
# Mapping unmapped IDs using the API... done.

# Log file path
log_file="3_logfile3.log"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$log_file"
}

# Input and output file paths
input_tsv="/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/chr_pos_to_ensp.tsv"
unique_ensp="/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/unique_ensp"
unmapped_file="/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/unmapped_ids.log"
output_json="/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/gffutils_ensp_to_uniprot_mappings.json"

batch_size=25  # Number of ENSP IDs per batch (adjustable)
failed_ids_dir="/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/failed_ids"
raw_dir="/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/raw"
successful_ids_dir="/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/successful_ids"
mkdir -p "$failed_ids_dir" "$raw_dir" "$successful_ids_dir"

# Extract unique ENSP IDs if the file doesn't already exist
if [ ! -f "$unique_ensp" ]; then
    log "Extracting unique ENSP IDs from $input_tsv"
    awk -F'\t' 'NR > 1 {print $6}' "$input_tsv" | sort -u > "$unique_ensp"
    log "Unique ENSP IDs written to $unique_ensp"
else
    log "File $unique_ensp already exists. Skipping extraction."
fi

"""
unique_ensp processed by 3_ensp_to_uniprot.py
Unmapped IDs written to $unmapped_file
"""

# Split unmapped ENSP IDs into batches
log "Splitting ENSP IDs into batches of size $batch_size"
split -l "$batch_size" "$unmapped_file" batch_

trap 'rm -f temp_results.json' EXIT

# Process each batch
batch_count=0
for batch_chunk in batch_*; do
    batch_count=$((batch_count + 1))

    # Log the batch being processed
    num_ids_in_batch=$(wc -l < "$batch_chunk")
    log "Processing batch $batch_count with $num_ids_in_batch ENSP IDs"

    # Create a comma-separated list of ENSP IDs for the batch
    batch_file=$(mktemp)
    log "Created temporary file $batch_file"
    paste -sd, "$batch_chunk" | tr -d '\r' | tr '\n' ',' | sed 's/,$//' > "$batch_file"

    # Submit the batch request to the UniProt API
    log "Submitting batch $batch_count to the API"
    response=$(curl --silent --request POST 'https://rest.uniprot.org/idmapping/run' \
        --form "ids=$(cat $batch_file)" \
        --form 'from="Ensembl_Protein"' \
        --form 'to="UniProtKB"')
    
    # Log the raw response from the API for debugging purposes
    log "Raw API response saved to debug_raw_results_batch_$batch_count.json"
    echo "$response" > "$raw_dir/debug_raw_results_batch_$batch_count".json

    # Extract jobId from the response
    jobId=$(echo "$response" | jq -r '.jobId')
    if [ -z "$jobId" ]; then
        log "Failed to retrieve jobId for batch: $(cat $batch_chunk)"
        continue
    else
        job_details=$(curl --silent "https://rest.uniprot.org/idmapping/details/$jobId")
        log "Job details for $jobId: $job_details"
    fi

    # Poll the job status
    log "Polling job status for jobId: $jobId"
    retry_count=0
    max_retries=3  # Retry limit for transient errors
    
    while true; do
        status_response=$(curl --silent "https://rest.uniprot.org/idmapping/status/$jobId")
        log "Job status response: $status_response"
        
        status=$(echo "$status_response" | jq -r '.jobStatus')
        if [ "$status" == "FINISHED" ]; then
            log "Job $jobId finished successfully."
            break
        elif [ "$status" == "FAILED" ]; then
            log "Job $jobId failed. Skipping batch."
            break
        elif [ "$status" == "ERROR" ]; then
            log "Job $jobId encountered an error: $(echo "$status_response" | jq -r '.errors[]?.message')"
            
            # Retry logic for transient errors
            retry_count=$((retry_count + 1))
            if [ "$retry_count" -le "$max_retries" ]; then
                log "Retrying job $jobId ($retry_count/$max_retries)..."
                sleep $((2 ** retry_count))  # Exponential backoff
            else
                log "Job $jobId failed after $max_retries retries. Marking batch as failed."
                echo "$batch_chunk" >> "$failed_ids_dir/failed_ids_batch_$batch_count.log"
                break
            fi
        else
            log "Job $jobId still in progress... waiting"
            sleep 5
        fi
    done

    # Fetch and parse results
    log "Fetching results for jobId: $jobId"
    result=$(curl --silent "https://rest.uniprot.org/idmapping/uniprotkb/results/$jobId")

    # Save failed IDs to a separate file
    failed_ids=$(echo "$result" | jq -r '.failedIds[]?' || true)
    if [ -n "$failed_ids" ]; then
        failed_ids_file="$failed_ids_dir/failed_ids_batch_$batch_count.log"
        echo "$failed_ids" > "$failed_ids_file"
        log "Failed IDs for batch $batch_count written to $failed_ids_file"
    else
        log "No failed IDs for batch $batch_count"
    fi

    if [ -z "$result" ]; then
        log "No results returned for batch: $(cat $batch_chunk)"
        continue
    fi

    # Process successful mappings
    successful_batch_file="$successful_ids_dir/successful_ids_batch_$batch_count.json"
    echo "$result" | jq -c '.results[]' > temp_results.json
    while read -r record; do
        ensp_id=$(echo "$record" | jq -r '.from')
        primaryAccession=$(echo "$record" | jq -r '.to.primaryAccession')
        if [ -n "$primaryAccession" ]; then
            log "Mapping found: $ensp_id -> $primaryAccession"
            echo "{\"from\": \"$ensp_id\", \"to\": \"$primaryAccession\"}" >> "$successful_batch_file"
            echo "{\"from\": \"$ensp_id\", \"to\": \"$primaryAccession\"}" >> "$output_json"
        fi
    done < temp_results.json

    log "Successful mappings for batch $batch_count written to $successful_batch_file"

    # Clean up temporary files
    rm -f "$batch_file" temp_results.json
done

: '
# Create an associative array to store mappings
declare -A mappings

# Ensure cleanup of temporary files on exit
trap 'rm -f "$temp_file" "$batch_file" batch_*' EXIT

# Process batches
for batch_chunk in batch_*; do
    log "Processing batch: $(basename "$batch_chunk")"
    paste -sd, "$batch_chunk" > "$batch_file"
    log "Submitting batch $(basename "$batch_chunk") with $(wc -l < "$batch_file") ENSP IDs"

    # Submit the batch request to the UniProt API
    response=$(curl --silent --request POST 'https://rest.uniprot.org/idmapping/run' \
        --form "ids=$(cat $batch_file)" \
        --form 'from="Ensembl_Protein"' \
        --form 'to="UniProtKB"')
    
    # Log the response from the API
    log "API response for batch $(basename "$batch_chunk"): $response"
    
    jobId=$(echo "$response" | jq -r '.jobId')
    if [ -z "$jobId" ]; then
        log "Failed to retrieve jobId for batch: $(cat $batch_chunk)"
        continue
    fi

    # Poll the job status
    log "Polling job status for jobId: $jobId"
    while true; do
        status_response=$(curl --silent "https://rest.uniprot.org/idmapping/status/$jobId")
        log "Job status response: $status_response"
        status=$(echo "$status_response" | jq -r '.jobStatus')
        
        if [ "$status" == "FINISHED" ]; then
            log "Job $jobId finished successfully."
            break
        elif [ "$status" == "FAILED" ]; then
            log "Job $jobId failed."
            continue 2
        else
            log "Job $jobId still in progress... waiting"
            sleep 5
        fi
    done

    # Fetch and parse results
    log "Fetching results for jobId: $jobId"
    result=$(curl --silent "https://rest.uniprot.org/idmapping/uniprotkb/results/$jobId")
    if [ -z "$result" ]; then
        log "No results returned for batch: $(cat $batch_chunk)"
        continue
    fi

    log "Processing results for batch: $(basename "$batch_chunk")"
    echo "$result" | jq -c '.results[]' > temp_results.json

    while read -r record; do
        log "Processing record: $record"
        ensp_id=$(echo "$record" | jq -r '.from')
        primaryAccession=$(echo "$record" | jq -r '.to.primaryAccession')
        
        if [ -z "$primaryAccession" ]; then
            log "No UniProt mapping found for ENSP ID: $ensp_id"
        else
            mappings["$ensp_id"]=$primaryAccession
        fi
    done < temp_results.json
done

# Write mappings to JSON
log "Writing mappings to $output_json"
echo "{" > "$output_json"
for key in "${!mappings[@]}"; do
    echo "  \"$key\": \"${mappings[$key]}\"," >> "$output_json"
done
sed -i '$ s/,$//' "$output_json"
echo "}" >> "$output_json"

log "Processing complete. Output saved to $output_json. Total mappings: ${#mappings[@]}"
'