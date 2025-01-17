'''
Input:
########
    * -i : A path to the CIVIC .csv file
    * -m : A path to the folder containing mapping files
    * -d : The name of the doid mapping file
    * -e : The name of the ensp to uniprot accession mapping file
    * -o : A path to the output folder

Output:
########
    * A .csv file with mutation data mapped to doid terms and uniprot accessions

Usage:
########
    * python map_civic_csv.py -h

    *Gives a description of the neccessary commands

    * python map_civic_csv.py -i <path/input_file.vcf> -m <path/mapping_folder> -d <doid_mapping_file_name> -e <enst_mapping_file_name> -o <path/>

    *Runs the script with the given input csv and outputs a csv with mutation mapped to doid terms and uniprot accession
'''

import argparse
from cmath import nan
import csv
import pandas as pd
from pathlib import Path
import re
import logging
import sys

sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent))
from utils import qc

logging.basicConfig(
    filename='map_civic.log',
    level=logging.DEBUG,
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    )

logging.info("Logging started--------------------------------")

def main(civic_csv, mapping_folder, doid_mapping_csv, enst_mapping_csv, output_folder):
    ##################################
    # Load the mapping files
    ##################################
    # Load in the TCGA mapping file to a mapping and cancer list
    doid_file_csv = mapping_folder + '/' + doid_mapping_csv
    #enst_file_csv = mapping_folder + '/' + enst_mapping_csv
    enst_file_csv = '/data/shared/repos/biomuta-old/downloads/glygen/' + enst_mapping_csv

    with open(doid_file_csv, "r") as doid_mapping_handle:
        doid_mapping = csv.reader(doid_mapping_handle)
        # Skip the header
        next(doid_mapping)

        # Set up the mapping dictionary
        doid_mapping_dict = {}
        # Populate the mapping dict
        for row in doid_mapping:
            doid_mapping_dict[row[0]] = row[1]
        
        for _,value in doid_mapping_dict.items():
            re.sub(r'NA', '', value)
    
    # Set up a list of cancers to iterate through
    cancer_list = []
    for _, value in doid_mapping_dict.items():
        if value not in cancer_list:
            cancer_list.append(value)
    
    # Load the ENSP to uniprot mapping file
    with open(enst_file_csv, "r") as enst_file_handle:
        enst_mapping = csv.reader(enst_file_handle, quoting=csv.QUOTE_ALL)
        # Skip the header
        next(enst_mapping)

        # Set up the mapping file dictionary
        ensp_mapping_dict = {}

        # Populate the mapping dictionary with keys as enst IDs and values as the gene symbol
        for row in enst_mapping:
            row[2] = re.sub(r'\.\d+$', '', row[2]) # Remove the period and digits from the transcript_id
            ensp_mapping_dict[row[2]] = row[1]
    
    ##################################
    # Load the civic csv file and map, then export
    ##################################
    civic_df = pd.read_csv(civic_csv, dtype=str)
    logging.info(f"Initial rows: {len(civic_df)}")


    # Map doid child to parent terms

    ## Strict equality check first
    civic_df['do_name'] = civic_df['CIViC Entity Disease'].map(doid_mapping_dict)

    ## Log the number of mapped items
    mapped_count = civic_df['do_name'].notna().sum()
    logging.info(f"Number of rows with mapped diseases after strict mapping: {mapped_count}")

    ## Then partial match check
    civic_df.loc[civic_df['do_name'].isna(), 'do_name'] = civic_df['CIViC Entity Disease'].apply(
        lambda x: map_partial_match(x, doid_mapping_dict)
        )
    mapped_count = civic_df['do_name'].notna().sum()
    logging.info(f"Number of rows with mapped diseases after partial match mapping: {mapped_count}")

    '''
    ## Log unmatched diseases
    for _, row in civic_df.iterrows():
        if pd.isna(row['do_name']):
            logging.warning(f"Unmatched disease: {row['CIViC Entity Disease']}")
    '''
            
    ## Convert 'NA' values to NaN
    civic_df['do_name'] = civic_df.apply(
        lambda row: convert_NA(row['do_name']),
        axis=1
        )
    logging.info(f"Rows converted to NA after DOID mapping: {civic_df['do_name'].isna().sum()}")


    '''
    # Check which entity diseases failed to map
    civic_diseases = civic_df['CIViC Entity Disease'].unique()
    doid_keys = list(doid_mapping_dict.keys())
    unmatched_diseases = [disease for disease in civic_diseases if disease not in doid_keys]
    logging.warning(f"Number of unmatched diseases: {len(unmatched_diseases)}")
    logging.warning(f"Unmatched diseases: {unmatched_diseases}")
    '''

    # Create a column that removes the dot notation from the ENST IDs in civic data
    civic_df['sample_name'] = ''
    civic_df['ENST'] = ''
    civic_df['ref_aa'] = ''
    civic_df['alt_aa'] = ''
    civic_df['aa_pos'] = ''
    civic_df['source'] = 'civic'
    civic_df['end_pos'] = ''

    # Check for indels and remove
    civic_df['ref_nt'] = civic_df['ref_nt'].apply(remove_indels)
    logging.info(f"Rows after removing indels: {len(civic_df)}")

    # Format the amino acid change and position
    logging.info('Formatting amino acid information')
    # amino acid changes to exclude
    logging.info(f"Rows before filtering on amino acid info: {len(civic_df)}")
    civic_df['amino_acid_info'] = civic_df['CIViC Variant Name'].apply(aa_format)
    logging.info(f"Rows with valid amino acid info: {civic_df['amino_acid_info'].notna().sum()}")
    civic_df.dropna(subset=['amino_acid_info'],inplace=True)
    civic_df[['ref_aa','alt_aa','aa_pos']] = pd.DataFrame(civic_df['amino_acid_info'].tolist(), index=civic_df.index)
    
    # Create a new column with only ENST ID to be used for mapping and separate AA notation
    civic_df['ENST'] = civic_df['Feature'].apply(lambda x: x.lstrip('_').split('.')[0])

    # Map ENST symbol to uniprot accession and store unmatched ENST in a string ready to be passed on to the UniProt API
    civic_df['uniprotkb_canonical_ac'] = civic_df['ENST'].map(ensp_mapping_dict)
    logging.info(f"Columns in dataframe: {civic_df.columns.tolist()}")
    logging.info(f"Rows with UniProt mapping: {civic_df['uniprotkb_canonical_ac'].notna().sum()}")
    logging.info(f"Rows that failed to map to UniProt: {civic_df['uniprotkb_canonical_ac'].isna().sum()}")
    unmatched_enst_str = ""  # Start with an empty string
    seen_enst = set()  # Set to track unique ENST IDs
    for _, row in civic_df.iterrows():
        enst_id = row['ENST']
        if pd.isna(row['uniprotkb_canonical_ac']) and qc.is_valid_enst(enst_id):
            if enst_id not in seen_enst:
                seen_enst.add(enst_id)  # Add to the set to track uniqueness
                if unmatched_enst_str:
                    unmatched_enst_str += ","
                unmatched_enst_str += enst_id  # Append the ENST ID to the string
    # Now unmatched_enst_str contains the unique comma-separated ENST IDs
    print(unmatched_enst_str)



    # Select and rename fields for integration with other sources
    final_fields = (
        'sample_name',
        'chr_id',
        'start_pos',
        'end_pos',
        'ref_nt',
        'alt_nt',
        'aa_pos',
        'ref_aa',
        'alt_aa',
        'do_name',
        'uniprotkb_canonical_ac',
        'source'
    )

    final_df = civic_df.loc[:, final_fields]

    #Final processing for the output df
    
    final_df['end_pos'] = final_df['start_pos']
    # Check which rows are being dropped
    missing_rows = final_df[final_df.isnull().any(axis=1)]
    missing_csv_path = '/data/shared/repos/biomuta-old/generated_datasets/civic/2025_01/missing.csv'
    missing_rows.to_csv(missing_csv_path, index = False)
    logging.info(f"Saved rows with NAs to file {missing_csv_path}")
    # Now drop rows with NA
    final_df.dropna(inplace=True)
    logging.info(f"After dropping NA: {len(final_df)}")
    final_df.drop_duplicates(keep='first',inplace=True)
    logging.info(f"After dropping duplicates: {len(final_df)}")

    mapped_new_file_path = output_folder + "/civic_missense_biomuta_v6.csv"
    logging.info(f"Final rows in output: {len(final_df)}")
    logging.info(f"Exporting mapped file to {mapped_new_file_path}")
    final_df.to_csv(mapped_new_file_path, index = False)

###############################
# Functions for formatting data
###############################

# Format the amino acid infomation
def aa_format(aa_info):
    # Define exceptions and additonal information to remove
    aa_exceptions = ['FRAMESHIFT','RS','rs','MUTATION','fs','FS','DEL','c.','DUP','HOM']
    aa_clean_up = ['BCR-ABL_','PML-RARA_','EM4-ALK_','ALK_Fusion_','HIP1-ALK_','FIP1L1-PDGFRA_','CD74-ROS1_','ETV6-NTRK3_']
    exception_flag = 0
    for info in aa_clean_up:
        aa_info = re.sub(info, '', str(aa_info))   
    for exception in aa_exceptions:
        if re.search(exception,aa_info):
            exception_flag = 1
            aa_list = [nan,nan,nan]

    # Format the amino acid change
    if exception_flag == 0:
        aa_list = re.findall(r'[A-Z\*]',aa_info)
        aa_position = re.findall(r'\d+',aa_info)
        aa_list.append(aa_position[0])

    # Account for additional outlier cases
    if len(aa_list) != 3:
        return [nan,nan,nan]
    else:
        return aa_list
    
def remove_indels(nt_info):
    if len(nt_info) > 1:
        nt_info = nan
    
    return nt_info

def convert_NA(NA_value):
    if NA_value in ['NA', 'None']:
        NA_value = nan
    
    return NA_value

def map_partial_match(value, mapping):
    # Convert 'NA' values to NaN using convert_NA(NA_value)
    value = convert_NA(value)
    
    # If value is NaN after conversion, just return it
    if pd.isna(value):
        return value

    # Normalize value by converting to lower case and removing underscores
    normalized_value = value.lower().replace('_', ' ')

    # Check for partial string match
    for key, mapped_value in mapping.items():
        # Normalize dictionary key
        normalized_key = key.lower().replace('_', ' ')
        if normalized_key in normalized_value:
            return mapped_value
    logging.debug(f"No match found for normalized_value: '{normalized_value}' in mapping keys.")
    return nan # Default if no match is found
    

                            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for civic mapping to doid and uniprot accessions.')
    parser.add_argument('--civic_csv', '-c',
                        help='An absolute path to the civic csv')
    parser.add_argument('--mapping_folder', '-m',
                        help='A path to the folder containing mapping files')                       
    parser.add_argument('--doid_mapping', '-d',
                        help='The name of the doid mapping file')
    parser.add_argument('--enst_mapping', '-e',
                        help='The name of the enst mapping file')
    parser.add_argument('--output_folder', '-o',
                        help='A path to the folder to export the mapped file')
    args = parser.parse_args()

    main(args.civic_csv, args.mapping_folder, args.doid_mapping, args.enst_mapping, args.output_folder)

# python3 map_civic_csv.py -c /data/shared/repos/biomuta-old/generated_datasets/civic/2025_01/civic_converted_mutations.csv -m /data/shared/repos/biomuta-old/pipeline/convert_step2/mapping -d civic_doid_mapping.csv -e human_protein_transcriptlocus.csv -o /data/shared/repos/biomuta-old/generated_datasets/civic/2025_01
