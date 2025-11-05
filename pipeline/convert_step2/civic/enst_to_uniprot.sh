#!/bin/bash

# Log file
log_file="enst_to_uniprot.log"
log() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$log_file" 2>&1 # Redirects stderr (file descriptor 2) to stdout (file descriptor 1). However, since the log is being appended to the file, this ensures that all log messages go to the log file, not to stdout.
}

# Get ENST IDs passed from map_civic_csv.py
enst_ids=$1
log "Received ENST string: $enst_ids"

# Submit ENST IDs to UniProt API
log "Submitting ENST IDs to the API"
response=$(curl --silent --request POST 'https://rest.uniprot.org/idmapping/run' \
    --form "ids=$enst_ids" \
    --form 'from="Ensembl_Transcript"' \
    --form 'to="UniProtKB"')

# Extract jobId from the response
jobId=$(echo "$response" | jq -r '.jobId')
if [ -z "$jobId" ]; then
    log "Failed to retrieve jobId"
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
        log "Job $jobId failed. Skipping."
        break
    elif [ "$status" == "ERROR" ]; then
        log "Job $jobId encountered an error: $(echo "$status_response" | jq -r '.errors[]?.message')"
        
        # Retry logic for transient errors
        retry_count=$((retry_count + 1))
        if [ "$retry_count" -le "$max_retries" ]; then
            log "Retrying job $jobId ($retry_count/$max_retries)..."
            sleep $((2 ** retry_count))  # Exponential backoff
        else
            log "Job $jobId failed after $max_retries retries."
            break
        fi
    else
        log "Job $jobId still in progress... waiting"
        sleep 5
    fi
done

# Fetch results
log "Fetching results for jobId: $jobId"
result=$(curl --silent "https://rest.uniprot.org/idmapping/uniprotkb/results/$jobId")

# Parse results and write to a JSON object in a loop to pass on back to Python
json_output="{}" # Initialize the JSON object as an empty object
log "Processing results"
while read -r record; do
    enst_id=$(echo "$record" | jq -r '.enst_id')
    primaryAccession=$(echo "$record" | jq -r '.primaryAccession')
    if [ -n "$primaryAccession" ]; then
        log "Mapping found: $enst_id -> $primaryAccession"
        # Update the JSON object
        json_output=$(jq --arg enst_id "$enst_id" \
                         --arg primaryAccession "$primaryAccession" \
                         '. + {($enst_id): $primaryAccession}' <<<"$json_output")
    else
        log "No primaryAccession found for $enst_id"
    fi
done < <(echo "$result" | jq -c '.results[] | {enst_id: .from, primaryAccession: .to.primaryAccession}')
echo "$json_output" # Output the JSON object to Python