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
    * A .csv file with mutation data and a .txt file with the .vcf headers

Usage:
########
    * python map_civic_csv.py -h

    *Gives a description of the neccessary commands

    * python map_civic_csv.py -i <path/input_file.vcf> -s <path/schema.json> -p output_prefix_name -o <path/>

    *Runs the script with the given input vcf and outputs a json file.

'''

import argparse
import csv
from distutils.log import INFO
import re
import json
import numpy as np
import pandas as pd


import os
import shutil
import sys


def main(civic_csv, mapping_folder, doid_mapping_csv, enst_mapping_csv, output_folder):
    ##################################
    # Load the mapping files
    ##################################
    # Load in the TCGA mapping file to a mapping and cancer list
    doid_file_csv = mapping_folder + '/' + doid_mapping_csv
    enst_file_csv = mapping_folder + '/' + enst_mapping_csv
    
    with open(doid_file_csv, "r") as doid_mapping_handle:
        doid_mapping = csv.reader(doid_mapping_handle)
        # Skip the header
        next(doid_mapping)

        # Set up the mapping dictionary
        doid_mapping_dict = {}

        # Populate the mapping dict
        for row in doid_mapping:
            doid_mapping_dict[row[0]] = row[1]
    
    # Set up a list of cancers to iterate through
    cancer_list = []
    for key, value in doid_mapping_dict.items():
        if value not in cancer_list:
            cancer_list.append(value)
    
    # Load the ENSP to uniprot mapping file.
    with open(enst_file_csv, "r") as enst_file_handle:
        enst_mapping = csv.reader(enst_file_handle, quoting=csv.QUOTE_ALL)
        # Skip the header.
        next(enst_mapping)

        # Set up the mapping file dictionary.
        ensp_mapping_dict = {}

        # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
        for row in enst_mapping:
            ensp_mapping_dict[row[2]] = row[0]
    
    ##################################
    # Load the civic csv file and map, then export
    ##################################
    civic_df = pd.read_csv(civic_csv, dtype=str)

    # Map doid child to parent terms
    civic_df['doid_name'] = civic_df['CIViC Entity Disease'].map(doid_mapping_dict)

    # Create a column that removes the dot notation from the ENST IDs in civic data
    civic_df['ENST'] = ''
    
    # Create a new column with only ENST ID to be used for mapping and separate AA notation
    for index, row in civic_df.iterrows():
        enst_info = str(row['Feature']).split('.')
        enst_id = enst_info[0]
        civic_df['ENST'] = enst_id

    # Map ENST symbol to uniprot accession
    civic_df['uniprotkb_canonical_ac'] = civic_df['ENST'].map(ensp_mapping_dict)

    # Select and rename fields for integration with other sources
    final_fields = [
        '#CHROM',
        'POS',
        'REF',
        'ALT',
        'CIViC Variant Name',
        'doid_name',
        'uniprotkb_canonical_ac'
    ]

    final_df = civic_df[final_fields]
    
    # Name the fields in the exported file according tot he BioMuta convention
    final_df.rename(columns={
        '#CHROM': 'chr_id', 
        'POS': 'chr_pos', 
        'REF': 'ref_nt', 
        'ALT': 'alt_nt', 
        'CIViC Variant Name': 'aa_change'
    }, inplace=True)

    print(final_df.columns)

    mapped_new_file_path = output_folder + "/mapped_tcga_mutations.csv"
    print("Exporting mapped file to " + mapped_new_file_path)
    final_df.to_csv(mapped_new_file_path, index = False)

                            

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

#python map_civic_csv.py -c /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/civic/civic_mutations_Mar_01_2022.csv -m /mnt/c/Users/caule/github_general/biomuta/downloader/mapping -d civic_doid_mapping.csv -e human_protein_transcriptlocus.csv -o .