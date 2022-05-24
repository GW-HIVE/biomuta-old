'''
Input:
########
    * -n : A path to the newly updated Biomuta mutation list
    * -p : A path to the previous update of Biomuta's mutation list
    * -o : A path to the output folder, where the data report will be exported

Output:
########
    * A data report comparing new AA sites to old AA sites for Biomuta

Usage:
########
    * python convert_icgc_vcf.py -h

    *Gives a description of the neccessary commands

    * python convert_icgc_vcf.py -i <path/input_file.vcf> -s <path/schema.json> -p output_prefix_name -o <path/>

    *Runs the script with the given input vcf and outputs a json file.

'''

import argparse
import csv
from distutils.log import INFO
import pandas as pd
import re
import json
import numpy as np



def main(tcga_mapping_file, uniprot_mapping_file, new_file):
    '''
    '''

    # Load in the TCGA mapping file to a mapping and cancer list
    with open(tcga_mapping_file, "r") as mapping_handle:
        mapping_csv = csv.reader(mapping_handle, delimiter="\t")
        # Skip the header
        next(mapping_csv)

        # Set up the mapping dictionary
        mapping_dict = {}

        # Populate the mapping dict
        for row in mapping_csv:
            doid_term = str(row[0]) + " / " + str(row[1])
            mapping_dict[doid_term] = row[2]

        print(mapping_dict)
    
    # Load the uniprot mapping file.
    with open(uniprot_mapping_file, "r") as uniprot_file_handle:
        uniprot_file_csv = csv.reader(uniprot_file_handle, quoting=csv.QUOTE_ALL)
        # Skip the header.
        next(uniprot_file_csv)

        # Set up the mapping file dictionary.
        uniprot_mapping_dict = {}

        # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
        for row in uniprot_file_csv:
            uniprot_mapping_dict[row[0]] = row[1]
    

    # Load the new file as a dataframe
    new_df = pd.read_csv("new_file")

    # Map TCGA study names to doid terms
    new_df['doid_term'] = new_df['project_short_name'].map(mapping_dict)

    # Map HUGO symbol to uniprot id
    new_df['uniprotkb_ac'] = new_df['Hugo_symbol'].map(uniprot_mapping_dict)

    

    


    #cancer_list = []
    #tcga_mapping_list = []
    #tcga_file = open(tcga_mapping_file, "r")
    #tcga_reader = csv.reader(tcga_file, delimiter="\t")
    #next(tcga_reader)
    #for row in tcga_reader:
        #tcga_mapping_list.append(row)
        # Assumes DOID is column 1 and cancer term is column 2 in the mapping file
        #cancer_list.append(str(row[0]) + " / " + str(row[1]))


    
    # Map the new file for cancer terms and uniprot IDs
    #new_mut_cancer_dict = {}
    #new_mut_gene_dict = {}
    #mut_list = []
    #new_mut_file = open(new_file, "w")
    #new_file_reader = csv.reader(new_mut_file)
    #next(new_file_reader)
    #for row in new_file_reader:
        # Assumes the TCGA project name is in the first column of the mut file and the 3 column of the mapping file
        #for mapping in tcga_mapping_list:
            #if row[1] == mapping[3]: 
                # Check if cancer type already exists in new mutation dict
                







    # Process each cancer separately
    for cancer in cancer_list:

        print("Processing " + cancer)

       # A set of lists to hold both new and old information
        #new_study_list = []
        #new_study_list_unique = []
        #old_study_list = []
        #old_study_list_unique = []
        #overlap_study_list = []

        # Get mutations per study in new update
        #with open(new_file) as nf:
            #for row in nf:
                #if re.search(study):
                    #new_study_list.append(row[new_file_pos_index] + row[new_file_aa_index])
        
        # Get mutations per study in old update
        #with open(old_file) as nf:
            #for row in nf:
                #if re.search(study):
                    #old_study_list.append(row[1] + row[2])

        # Compare list od AA per study
        #for aa in new_study_list:
            






        # Get counts in old update


                            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for civic vcf to csv convertor.')
    parser.add_argument('--tcga_mapping_file', '-m',
                        help='An absolute path to the TCGA study to cancer mapping file')
    args = parser.parse_args()

    main(args.tcga_mapping_file)