import csv
import logging
import os

# Configure logging
log_dir = "/home/maria.kim/logs/2_convert"
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, "ensp_mapping_process.log")

# Set up logging configuration
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(log_file),
        logging.StreamHandler()  # This keeps console output as well
    ]
)

logger = logging.getLogger(__name__)

# File paths
input_file = "/data/shared/repos/biomuta-old/generated_datasets/current/mapping_ids/chr_pos_to_ensp.tsv"
mapping_file = "/data/shared/repos/biomuta-old/downloads/glygen/human_protein_transcriptlocus.csv"
output_file = "/data/shared/repos/biomuta-old/generated_datasets/current/mapping_ids/ensp_to_uniprot_from_glygen.json"
unmapped_file = "/data/shared/repos/biomuta-old/generated_datasets/current/mapping_ids/unmapped_ids_by_glygen.log"

def read_input_ids(input_path):
    """
    Generator function that extracts ENSP IDs from TSV file.
    Reads the 'ensp' column from the TSV file.
    """
    with open(input_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        processed_count = 0
        for row in reader:
            processed_count += 1
            if processed_count % 10000 == 0:
                logger.info(f"Processed {processed_count} rows from input TSV...")
            ensp_id = row["ensp"].strip()
            if ensp_id:  # Only yield non-empty ENSP IDs
                yield ensp_id

def process_mapping_file(mapping_path, ensp_set):
    """
    Generator function that processes mapping file and yields ENSP-UniProt pairs.
    """
    with open(mapping_path, 'r') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i % 10000 == 0:
                logger.info(f"Processed {i} rows from mapping file...")
            peptide_id = row["peptide_id"].split('.')[0]
            uniprot_ac = row["uniprotkb_canonical_ac"]
            if peptide_id in ensp_set:
                yield peptide_id, uniprot_ac

def write_output(output_path, mapping_generator):
    """
    Write ENSP-UniProt mappings incrementally to JSON.
    """
    logger.info("Writing output to JSON file...")
    with open(output_path, 'w') as f:
        f.write("{\n")
        first = True
        count = 0
        for peptide_id, uniprot_ac in mapping_generator:
            if not first:
                f.write(",\n")
            f.write(f'    "{peptide_id}": "{uniprot_ac}"')
            first = False
            count += 1
            if count % 10000 == 0:
                logger.info(f"Written {count} mappings to JSON...")
        f.write("\n}\n")
    logger.info(f"Finished writing {count} mappings to JSON.")

def log_unmapped_ids(input_ids, mapped_ids, log_path):
    """Log unmapped ENSP IDs to a file."""
    unmapped_ids = input_ids - mapped_ids
    logger.info(f"Logging {len(unmapped_ids)} unmapped IDs...")
    with open(log_path, 'w') as f:
        for unmapped_id in unmapped_ids:
            f.write(f"{unmapped_id}\n")
    logger.info("Unmapped IDs logging completed.")

# Main execution
if __name__ == "__main__":
    logger.info("Starting ENSP to UniProt mapping process...")
    
    logger.info("Reading input ENSP IDs from TSV file...")
    ensp_ids = set(read_input_ids(input_file))
    logger.info(f"Loaded {len(ensp_ids)} unique ENSP IDs from input TSV file.")

    logger.info("Processing mapping file...")
    mapped_ids = set()

    mapping_gen = ((peptide_id, uniprot_ac) for peptide_id, uniprot_ac in process_mapping_file(mapping_file, ensp_ids))
    mapping_gen_for_logging = ((peptide_id, uniprot_ac) for peptide_id, uniprot_ac in process_mapping_file(mapping_file, ensp_ids))

    for peptide_id, uniprot_ac in mapping_gen_for_logging:
        mapped_ids.add(peptide_id)

    write_output(output_file, mapping_gen)

    log_unmapped_ids(ensp_ids, mapped_ids, unmapped_file)

    logger.info(f"Mapping completed. Output written to {output_file}")
    logger.info(f"Unmapped IDs logged to {unmapped_file}")
    logger.info(f"Process log saved to {log_file}")
    logger.info("Process completed successfully!")