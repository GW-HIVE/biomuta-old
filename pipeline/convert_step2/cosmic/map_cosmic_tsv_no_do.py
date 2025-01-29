'''
Input:
########
    * -i : A path to the cosmic tsv mutation file
    * -m : A path to the folder containing mapping files
    * -e : The name of the enst to uniprot accession mapping file
    * -o : A path to the the folder to export the final mapped mutations


Output:
########
    * A mutation file with COSMIC mutations mapped to doid terms and uniprot accessions

Usage:
########
    * map_cosmic_tsv -h

    *Gives a description of the neccessary commands

    * python map_cosmic_tsv.py -i <path/cosmic_file_name.tsv> -m <path/mapping_folder> -e <enst_mapping_file_name> -o <path/output_folder>

    *Runs the script with the given input tsv and outputs a csv with Biomuta formatting.

'''

import argparse
from cmath import nan
import csv
import logging
import pandas as pd
import re

logging.basicConfig(filename="cosmic.log",
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    filemode='w',
                    level=logging.INFO)
logging.info("Logging started ------------------------")

invalid_aa_count = 0 # Count the number of 'p.?' in the column 'AA_MUT_SYNTAX'
invalid_entries = set() # Keep track of gene names missing their ENST ID
missing_enst_count = 0

def main(cosmic_tsv, mapping_folder, enst_mapping_csv, output_folder):
    global invalid_aa_count, invalid_entries, missing_enst_count
    ##################################
    # Load the mapping files
    ##################################
    enst_file_csv = mapping_folder + '/' + enst_mapping_csv
    
    # Load the ENST to uniprot mapping file.
    logging.info("Loading ENST mapping file...")
    with open(enst_file_csv, "r") as enst_file_handle:
        enst_mapping = csv.reader(enst_file_handle, quoting=csv.QUOTE_ALL)
        # Skip the header.
        next(enst_mapping)

        # Set up the mapping file dictionary.
        ensp_mapping_dict = {}

        # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
        for row in enst_mapping:
            ensp_mapping_dict[row[2]] = row[1]
    
    ##################################
    # Load the cosmic tsv file and map, then export
    ##################################
    columns_from_cosmic = [
        'CHROMOSOME',
        'GENOMIC_MUT_START', 
        'GENOMIC_MUT_STOP', 
        'GENOMIC_WT_ALLELE_SEQ', 
        'GENOMIC_MUT_ALLELE_SEQ', 
        'GENE_NAME',
        'AA_MUT_SYNTAX',
        ]
    
    logging.info(f"Loading mutations from COSMIC tsv mutation file for fields: {columns_from_cosmic}")

    cosmic_df_iterator = pd.read_csv(cosmic_tsv, usecols=columns_from_cosmic, dtype=str, sep='\t',chunksize=1000000)

    for i, cosmic_df in enumerate(cosmic_df_iterator):

        logging.info(f"Starting chunk {i}")
    
        # Create new fields for reformatted data: ENST, genome location, AA mutation, nucleotide mutation
        cosmic_df['ENST'] = ''
        cosmic_df.rename(columns = {'CHROMOSOME':'chr_id'},inplace=True)
        cosmic_df.rename(columns = {'GENOMIC_MUT_START':'start_pos'},inplace=True)
        cosmic_df.rename(columns = {'GENOMIC_MUT_STOP':'end_pos'},inplace=True)
        cosmic_df['ref_aa'] = ''
        cosmic_df['alt_aa'] = ''
        cosmic_df['aa_pos'] = ''
        cosmic_df.rename(columns = {'GENOMIC_WT_ALLELE_SEQ':'ref_nt'},inplace=True)
        cosmic_df.rename(columns = {'GENOMIC_MUT_ALLELE_SEQ':'alt_nt'},inplace=True)
    
        # Drop rows with missing information
        cosmic_df = drop_na_by_column(cosmic_df)
    
        # Create a new column with only ENST ID to be used for mapping. Also separate the AA notation, genome locations, and nucleotide change
    
        # Format the amino acid change and position
        logging.info('Formatting amino acid information')
        cosmic_df = cosmic_df.apply(aa_format, axis=1)
        logging.info(f"Cumulative count of 'p.?': {invalid_aa_count}")
        cosmic_df.dropna(subset=['amino_acid_info'],inplace=True)
        cosmic_df[['ref_aa','alt_aa','aa_pos']] = pd.DataFrame(cosmic_df['amino_acid_info'].tolist(), index=cosmic_df.index)
        
        # Map ENST symbol to uniprot accession
        logging.info('Mapping ENST IDs to uniprot accession')
        ## Apply the extract_enst function with the full row passed as an argument
        cosmic_df = cosmic_df.apply(extract_enst, axis=1)
        # Log the cumulative number of invalid entries and the cumulative number of unique genes
        logging.info(f"Cumulative number of invalid 'GENE_NAME' entries: {missing_enst_count}")
        logging.info(f"Cumulative number of unique 'GENE_NAME' entries: {len(invalid_entries)}")
        cosmic_df['uniprotkb_canonical_ac'] = cosmic_df['ENST'].map(ensp_mapping_dict)
    
        cosmic_df.dropna(subset=['ref_nt', 'ref_aa', 'chr_id'],inplace=True)
     
        # Select and rename fields for integration with other sources
        final_fields = (
            'chr_id',
            'start_pos',
            'end_pos',
            'ref_nt',
            'alt_nt',
            'aa_pos',
            'ref_aa',
            'alt_aa',
            'uniprotkb_canonical_ac'
        )
    
        final_df = cosmic_df.loc[:, final_fields]
        
        # Remove rows that did not map to uniprot canonical transcripts and duplicates
        final_df.dropna(subset=['uniprotkb_canonical_ac'],inplace=True)
        final_df.drop_duplicates(keep='first',inplace=True)
        
        # How to handle the df chunks to process
        mode = 'w' if i == 0 else 'a'
        header = i == 0
        

        mapped_new_file_path = output_folder + "/cosmic_missense_biomuta_v5.csv"
        logging.info(f"Adding processed data to {mapped_new_file_path}")
        final_df.to_csv(mapped_new_file_path, index = False, header=header, mode=mode)

        logging.info(f"Chunk number {i} completed")

    # Log final totals after all chunks are processed
    logging.info(f"Final total count of 'p.?': {invalid_aa_count}")
    logging.info(f"Final total number of invalid 'GENE_NAME' entries: {missing_enst_count}")
    logging.info(f"Final total number of unique 'GENE_NAME' entries: {len(invalid_entries)}")
    logging.info(f"Total number of invalid rows: {cosmic_df['is_valid'].value_counts()}")


###############################
# Functions for formatting data
###############################

# Drop NA column by column
import pandas as pd

def drop_na_by_column(df):
    """
    Drops rows with NA values one column at a time and logs the remaining entries.
    
    Args:
        df (pd.DataFrame): The DataFrame to process.
        
    Returns:
        pd.DataFrame: The modified DataFrame with rows dropped for each column.
    """
    for col in df.columns:
        initial_count = len(df)
        df.dropna(subset=[col], inplace=True)
        remaining_count = len(df)
        logging.info(f"After dropping rows with NA in '{col}': {remaining_count} entries remain (dropped {initial_count - remaining_count} rows).")
    return df

# Format the amino acid infomation
def aa_format(row):
    global invalid_aa_count
    if row['AA_MUT_SYNTAX'] == "p.?":
        invalid_aa_count += 1
        row['amino_acid_info'] = [nan, nan, nan]
        row['is_valid'] = False
        return row
    else:
        aa_list = re.findall(r'[A-Z\*]',row['AA_MUT_SYNTAX'])
        aa_position = re.findall(r'\d+',row['AA_MUT_SYNTAX'])
        aa_list.append(aa_position[0])
        if len(aa_list) != 3:
            row['amino_acid_info'] = [nan, nan, nan]
            row['is_valid'] = False
            return row
        else:
            row['is_valid'] = True
            return row
        
def extract_enst(row):
    global invalid_entries, missing_enst_count
    parts = row['GENE_NAME'].split('_')
    if len(parts) > 1: # Check if there are at least two parts after splitting by '_'
        row['ENST'] = parts[1]
        return row
    else: # Track the unique invalid gene names (missing ENST)
        invalid_entries.add(row['GENE_NAME'])
        missing_enst_count += 1
        row['is_valid'] = False
        return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for civic mapping to doid and uniprot accessions.')
    parser.add_argument('--cosmic_tsv', '-i',
                        help='An absolute path to the cosmic tsv')
    parser.add_argument('--mapping_folder', '-m',
                        help='A path to the folder containing mapping files')                       
    parser.add_argument('--enst_mapping', '-e',
                        help='The name of the enst mapping file')
    parser.add_argument('--output_folder', '-o',
                        help='A path to the folder to export the mapped file')
    args = parser.parse_args()

    main(args.cosmic_tsv, args.mapping_folder, args.enst_mapping, args.output_folder)
        
#python3 map_cosmic_tsv_no_do.py -i /data/shared/repos/biomuta-old/downloads/cosmic/cosmic_snps_narrow.tsv -m /data/shared/repos/biomuta-old/pipeline/convert_step2/mapping -e human_protein_transcriptlocus.csv -o /data/shared/biomuta/generated/datasets/cosmic/2025_01
