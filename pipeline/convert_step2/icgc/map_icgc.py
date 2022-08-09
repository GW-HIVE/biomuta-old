'''
Input:
########
    * -i : A path to the ICGC .csv file
    * -m : A path to the folder containing mapping files
    * -d : The name of the doid mapping file
    * -e : The name of the ensp to uniprot accession mapping file
    * -o : A path to the output folder

Output:
########
    * A .csv file with mutation data formatted to the biomuta field structure

Usage:
########
    * python map_icgc.py -h

    *Gives a description of the neccessary commands

    * python map_icgc.py -i <path/input_file.vcf> -m <path/> -d doid_mapping_file.csv -e enst_mapping_file.csv -o <path/>

    *Runs the script with the given csv file and outputs a csv file formatted for the final biomuta master file

'''

import argparse
from cmath import nan
import csv
import re
import pandas as pd
import numpy as np

def main(icgc_csv, mapping_folder, doid_mapping_csv, enst_mapping_csv, output_folder):
    ##################################
    # Load the mapping files
    ##################################
    # Load in the TCGA mapping file to a mapping and cancer list
    doid_file_csv = mapping_folder + '/' + doid_mapping_csv
    enst_file_csv = mapping_folder + '/' + enst_mapping_csv
    
    with open(doid_file_csv, "r") as doid_mapping_handle:
        doid_mapping = csv.reader(doid_mapping_handle, delimiter='\t')
        # Skip the header
        next(doid_mapping)

        # Set up the mapping dictionary
        doid_mapping_dict = {}

        # Populate the mapping dict with the TCGA study name as keys and the doid name as values
        for row in doid_mapping:
            tcga_study = re.sub('TCGA-', '', row[2])
            doid_term = str(row[0]) + " / " + str(row[1])
            doid_mapping_dict[tcga_study] = doid_term
    
    # Load the ENSP to uniprot mapping file.
    with open(enst_file_csv, "r") as enst_file_handle:
        enst_mapping = csv.reader(enst_file_handle, quoting=csv.QUOTE_ALL)
        # Skip the header.
        next(enst_mapping)

        # Set up the mapping file dictionary.
        ensp_mapping_dict = {}

        # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
        for row in enst_mapping:
            ensp_mapping_dict[row[2]] = row[1]
    
    columns_from_icgc = [
        'chr_id',
        'start_pos', 
        'ref_nt', 
        'alt_nt', 
        'aa_mutation', 
        'project_code', 
        'transcript_affected'
        ]
    
    print("Loading mutations from ICGC csv mutation file for fields: " + str(columns_from_icgc))

    icgc_df_iterator = pd.read_csv(icgc_csv, usecols=columns_from_icgc, dtype=str, chunksize=1000000)

    for i, icgc_df in enumerate(icgc_df_iterator):

        print("Starting chunk " + str(i))
    
        # Create new fields for reformatted data: ENST, genome location, AA mutation, nuceotide mutation
        icgc_df['ref_aa'] = ''
        icgc_df['alt_aa'] = ''
        icgc_df['aa_pos'] = ''
        icgc_df['do_name'] = ''
        icgc_df['tcga_study'] = ''
        icgc_df['end_pos'] = ''

        # Map dTCGA study terms to doid terms
        icgc_df['tcga_study'] = icgc_df['project_code'].apply(clean_tcga_code)
        icgc_df['do_name'] = icgc_df['tcga_study'].map(doid_mapping_dict)
        icgc_df.dropna(subset=['do_name'],inplace=True)

        # Format the amino acid change and position
        print('Formatting amino acid information')
        # Remove entries with no aa change
        icgc_df.dropna(subset=['aa_mutation'],inplace=True)
        icgc_df['amino_acid_info'] = icgc_df['aa_mutation'].apply(aa_format)
        # Remove entries with atypical aa change
        icgc_df.dropna(subset=['amino_acid_info'],inplace=True)
        icgc_df[['ref_aa','alt_aa','aa_pos']] = pd.DataFrame(icgc_df['amino_acid_info'].tolist(), index=icgc_df.index)
    
        # Format the genomic location
        icgc_df['end_pos'] = icgc_df['start_pos']
        
        # Map ENST symbol to uniprot accession
        print('Mapping ENST IDs to uniprot accession')
        icgc_df['uniprotkb_canonical_ac'] = icgc_df['transcript_affected'].map(ensp_mapping_dict)
        icgc_df.dropna(subset=['uniprotkb_canonical_ac'],inplace=True)

        # Select and rename fields for integration with other sources
        final_fields_tuple = (
            'chr_id',
            'start_pos',
            'end_pos',
            'ref_nt',
            'alt_nt',
            'aa_pos',
            'ref_aa',
            'alt_aa',
            'do_name',
            'uniprotkb_canonical_ac'
        )
    
        final_df = icgc_df.loc[:, final_fields_tuple]
    
        # Final processing for the output df
        final_df['sample_name'] = ''
        final_df['source'] = 'icgc'
        final_df.drop_duplicates(keep='first',inplace=True)
        
        # How to handle the df chunks to process
        mode = 'w' if i == 0 else 'a'
        header = i == 0
        

        mapped_new_file_path = output_folder + "/icgc_missense_biomuta_v5.csv"
        print("Adding processed data to " + mapped_new_file_path)
        final_df.to_csv(mapped_new_file_path, index = False, header=header, mode=mode)

        print("Chunk number " + str(i) + " completed")    

###############################
# Functions for formatting data
###############################

# Format the amino acid infomation
def aa_format(aa_info):
    # Grab all aa information
    aa_list = re.findall(r'[A-Z\*]',aa_info)
    aa_position = re.findall(r'\d+',aa_info)
    aa_list.append(aa_position[0])

    # Account for additional outlier cases
    if len(aa_list) != 3:
        return [nan,nan,nan]
    else:
        return aa_list

def clean_tcga_code(project_code):
    tcga_code = re.sub(r'OCCURRENCE=', '', project_code)
    tcga_code = re.sub(r'\-.*', '', tcga_code)
    return tcga_code

def convert_NA(NA_value):
    if NA_value == 'NA':
        NA_value = nan
    
    return NA_value
    

                            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for icgc mapping to doid and uniprot accessions.')
    parser.add_argument('--civic_csv', '-c',
                        help='An absolute path to the icgc csv')
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

#python map_icgc.py -c /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/icgc/icgc_converted_mutations.csv -m /mnt/c/Users/caule/github_general/biomuta/pipeline/convert_step2/mapping -d tcga_doid_mapping.csv -e human_protein_transcriptlocus.csv -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled 