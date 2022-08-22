'''
Input:
########
    * -n : A path to the newly updated Biomuta mutation list


Output:
########
    * A data report comparing new AA sites to old AA sites for Biomuta

Usage:
########
    * python process_tcga_download.py -h

    *Gives a description of the neccessary commands

    * python process_tcga_download.py -i <path/input_file.vcf> -s <path/schema.json> -p output_prefix_name -o <path/>

    *Runs the script with the given input vcf and outputs a json file.

'''

import argparse
import csv
from distutils.log import INFO
from genericpath import exists
import pandas as pd
import re
import json
import numpy as np



def main(mapping_folder, cosmic_cancer_types, doid_mapping):
    
    doid_mapping_path = mapping_folder + '/' + doid_mapping
    cosmic_cancers_path = mapping_folder + '/' + cosmic_cancer_types

    doid_handle = open(doid_mapping_path, "r")
    doid_cancers = csv.reader(doid_handle, delimiter="\t")

    cosmic_handle = open(cosmic_cancers_path, "r")
    cosmic_cancers = csv.reader(cosmic_handle)
    
    # What fields are being mapped from each input file?
    cosmic_mapping_field_index = 0

    # The field with terms that have been mapped to DOID
    doid_matched_terms_index = 11

    # The field with doid paretn terms
    doid_mapping_field_index = 3
        
    # grab headers from input to create output headers
    cosmic_header = next(cosmic_cancers)
    cosmic_field_name = cosmic_header[cosmic_mapping_field_index]
    doid_header = next(doid_cancers)
    doid_field_name = doid_header[doid_mapping_field_index]

    print('Mapping Cosmic: ' + cosmic_field_name + ' to DOID: ' + doid_field_name)
    
    out_header = [cosmic_field_name,doid_field_name]

    # Create the doid mapping information dictionary
    doid_dict = {}
    for doid_row in doid_cancers:
        doid_dict[doid_row[doid_matched_terms_index]] = doid_row[doid_mapping_field_index]




    mapping_dict = {}
    
    # Create the mapping dictionary with cosmic cancer type as keys and doid terms as values
    for cosmic_row in cosmic_cancers:

        if cosmic_row[cosmic_mapping_field_index] == 'NS':
            mapping_dict['NS'] = ['']
            continue

        print(cosmic_row)
        
        mapping_logged = False

        for key, value in doid_dict.items():
            
            # When the mapping exists, go to next cancer type
            if mapping_logged == True:
                break
            
            if re.search(str(cosmic_row[cosmic_mapping_field_index]), key):
                print('Mapping ' + cosmic_row[cosmic_mapping_field_index] + ' to ' + value)
                mapping_logged = True
                mapping_dict[cosmic_row[cosmic_mapping_field_index]] = value
    
    print(mapping_dict)

    cosmic_doid_mapping = mapping_folder + '/cosmic_doid_mapping.csv'

    with open(cosmic_doid_mapping, 'w') as out_file:
        writer = csv.DictWriter(out_file, fieldnames=out_header)
        writer.writeheader()
        for mapping in mapping_dict:
            writer.writerow(mapping)
            
        

    





                            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for creating the cosmic cancer type to doid mapping file')
    parser.add_argument('--mapping_folder', '-m',
                        help='An absolute path to the folder containing mapping files')
    parser.add_argument('--cosmic_cancer_types', '-c',
                        help='An absolute path to the new cosmic cancer types file')
    parser.add_argument('--doid_mapping', '-d',
                        help='An absolute path to the doid mapping file')
    args = parser.parse_args()

    main(args.mapping_folder, args.cosmic_cancer_types, args.doid_mapping)


    # Example run:

#python make_cosmic_doid_cancer_mapping.py -c -d