#!/bin/bash

# API only processes 25 IDs at a time for some reason
# This script processes unmapped ENSP IDs from the previous GlyGen mapping step
# and maps them to UniProt IDs using the UniProt API

# Log file path
log_file="/home/maria.kim/logs/2_convert/3_2_ensp_to_uniprot_api_updated.log"

log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$log_file"
}

# Input and output file paths
unmapped_file="/data/shared/repos/biomuta-old/generated_datasets/current/mapping_ids/unmapped_ids_by_glygen.log"
output_json="/data/shared/repos/biomuta-old/generated_datasets/current/mapping_ids/ensp_to_uniprot_from_api.json"

batch_size=25  # Number of ENSP IDs per batch (adjustable)
failed_ids_dir="/data/shared/repos/biomuta-old/generated_datasets/current/mapping_ids/failed_ids"
raw_dir="/data/shared/repos/biomuta-old/generated_datasets/current/mapping_ids/raw"
successful_ids_dir="/data/shared/repos/biomuta-old/generated_datasets/current/mapping_ids/successful_ids"
mkdir -p "$failed_ids_dir" "$raw_dir" "$successful_ids_dir"

# Check if unmapped file exists
if [ ! -f "$unmapped_file" ]; then
    log "Error: Unmapped IDs file $unmapped_file not found. Please run the previous GlyGen mapping step first."
    exit 1
fi

# Split unmapped ENSP IDs into batches
log "Splitting ENSP IDs into batches of size $batch_size"
split -l "$batch_size" "$unmapped_file" batch_

trap 'rm -f temp_results.json; rm -f batch_*' EXIT

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
    poll_count=0
    max_polls=240  # Maximum 20 minutes of polling (240 * 5 seconds)
    
    while true; do
        status_response=$(curl --silent "https://rest.uniprot.org/idmapping/status/$jobId")
        log "Job status response: $status_response"
        
        status=$(echo "$status_response" | jq -r '.jobStatus')
        poll_count=$((poll_count + 1))
        
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
                cat "$batch_chunk" >> "$failed_ids_dir/failed_ids_batch_$batch_count.log"
                break
            fi
        elif [ "$poll_count" -ge "$max_polls" ]; then
            log "Job $jobId has been stuck in '$status' status for too long (40+ minutes). Marking batch as failed and moving on."
            cat "$batch_chunk" >> "$failed_ids_dir/failed_ids_batch_$batch_count.log"
            break
        else
            log "Job $jobId still in progress (status: $status)... waiting (poll $poll_count/$max_polls)"
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

    # Initialize a temporary file for the dictionary if it doesn't exist
    if [ ! -f "$output_json" ]; then
        echo "{}" > "$output_json"
    fi

    successful_mappings=$(mktemp)
    echo "{}" > "$successful_mappings" # Temporary file to accumulate results

    # Process successful mappings
    log "Processing successful mappings for batch $batch_count"
    successful_batch_file="$successful_ids_dir/successful_ids_batch_$batch_count.json"
    echo "$result" | jq -c '.results[]' > temp_results.json
    
    # Transform and accumulate key-value pairs
    while read -r record; do
        ensp_id=$(echo "$record" | jq -r '.from')
        primaryAccession=$(echo "$record" | jq -r '.to.primaryAccession')
        if [ -n "$primaryAccession" ]; then
            log "Mapping found: $ensp_id -> $primaryAccession"
            jq --arg key "$ensp_id" --arg value "$primaryAccession" \ '. + {($key): $value}' "$successful_mappings" > tmp.json && mv tmp.json "$successful_mappings"
        fi
    done < temp_results.json

    # Merge the batch results with the main output file
    jq -s '.[0] * .[1]' "$output_json" "$successful_mappings" > tmp.json && mv tmp.json "$output_json"
    log "Merged successful mappings for batch $batch_count into $output_json"

    # Clean up temporary files
    rm -f "$batch_chunk" "$batch_file" temp_results.json "$successful_mappings"
done

log "Script completed. Results written to $output_json"