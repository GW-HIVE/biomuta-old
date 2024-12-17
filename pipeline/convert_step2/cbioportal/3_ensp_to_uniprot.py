import csv
import json

# File paths
input_file = "/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/unique_ensp"
mapping_file = "/data/shared/repos/biomuta-old/downloads/glygen/human_protein_transcriptlocus.csv"
output_file = "/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/ensp_to_uniprot.json"
unmapped_file = "/data/shared/repos/biomuta-old/generated_datasets/2024_10_22/mapping_ids/unmapped_ids.log"

def read_input_ids(input_path):
    """Generator to yield ENSP IDs from input file."""
    with open(input_path, 'r') as f:
        for line in f:
            yield line.strip()

def process_mapping_file(mapping_path, ensp_set):
    """Generator to process mapping file and yield ENSP-UniProt pairs."""
    with open(mapping_path, 'r') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i % 10000 == 0:
                print(f"Processed {i} rows from mapping file...")
            peptide_id = row["peptide_id"].split('.')[0]
            uniprot_ac = row["uniprotkb_canonical_ac"].split('-')[0]
            if peptide_id in ensp_set:
                yield peptide_id, uniprot_ac

def write_output(output_path, mapping_generator):
    """Write ENSP-UniProt mappings incrementally to JSON."""
    print("Writing output to JSON file...")
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
                print(f"Written {count} mappings to JSON...")
        f.write("\n}\n")
    print(f"Finished writing {count} mappings to JSON.")

def log_unmapped_ids(input_ids, mapped_ids, log_path):
    """Log unmapped ENSP IDs to a file."""
    unmapped_ids = input_ids - mapped_ids
    print(f"Logging {len(unmapped_ids)} unmapped IDs...")
    with open(log_path, 'w') as f:
        for unmapped_id in unmapped_ids:
            f.write(f"{unmapped_id}\n")
    print("Unmapped IDs logging completed.")

# Main execution
if __name__ == "__main__":
    print("Reading input ENSP IDs...")
    ensp_ids = set(read_input_ids(input_file))
    print(f"Loaded {len(ensp_ids)} ENSP IDs from input file.")

    print("Processing mapping file...")
    mapped_ids = set()

    mapping_gen = ((peptide_id, uniprot_ac) for peptide_id, uniprot_ac in process_mapping_file(mapping_file, ensp_ids))
    mapping_gen_for_logging = ((peptide_id, uniprot_ac) for peptide_id, uniprot_ac in process_mapping_file(mapping_file, ensp_ids))

    for peptide_id, uniprot_ac in mapping_gen_for_logging:
        mapped_ids.add(peptide_id)

    write_output(output_file, mapping_gen)

    log_unmapped_ids(ensp_ids, mapped_ids, unmapped_file)

    print(f"Mapping completed. Output written to {output_file}")
    print(f"Unmapped IDs logged to {unmapped_file}")
