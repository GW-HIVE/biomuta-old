'''
Input:
########
    * -i : A path to the CIVIC .vcf file
    * -p : A prefix used for naming the output files
    * -o : A path to the output folder, where the data csv and vcf headers will go

Output:
########
    * A .csv file with mutation data and a .txt file with the .vcf headers

Usage:
########
    * python convert_civic_vcf.py -h

    *Gives a description of the neccessary commands

    * python convert_civic_vcf -i <path/input_file.vcf> -s <path/schema.json> -p output_prefix_name -o <path/>

    *Runs the script with the given input vcf and outputs a json file.

'''

import argparse
import csv
import re
import json
import numpy as np


import os
import shutil
import sys


def main(input_vcf_file, schema, output_prefix, output_folder):
    '''
    Loads an input .vcf and converts rows to csv format 
    '''
    
    # A list to hold the vcf headers
    vcf_headers = []

    # A list to hold the mutation data
    data = []    

    # Separate the VCF headers and mutation data
    with open(input_vcf_file, 'r') as vcf:
        
        reader = csv.reader(vcf, delimiter="\t")
        
        # Go through the VCF by row
        for row in reader:

            # Pull out VCF headers
            if re.search(r'##', str(row)):
                vcf_headers.append(row)
                
                # Record number of subfields in the INFO column
                if re.search(r'##INFO=<ID=CSQ', str(row)):
                    sub_field_total = str(row).count('|')
                    #print("The INFO column contains " + str(sub_field_total) + " fields")


            
            
            # Pull out header for mutation data
            elif re.search(r'#CHROM', str(row)):
                
                mut_fields = row

            #elif re.search(r'fileformat', str(row)):
                
            #    print(row)
            
            else:
                data.append(row)
        
    # Start a counter and end count used for troubleshooting
    count = 0 
    stop_count = 3
    
    # A list of mutation dictionaries that will be converted into the output json
    dict_list = []
    
    # Convert each line in the mutation data to json format
    for line in data: 

        # Count the number of sub_fields in the info column
        sub_field_count = str(line).count('|')

        #print("Line " + str(count) + " contains " + str(int(sub_field_count/sub_field_total)) + " annotations")
        #print(str(sub_field_count/sub_field_total))

        # A dictionary to hold mutation data
        mut_dict = {}

        # A list and dictionary to hold the info subfields
        annotation_count = int(sub_field_count/sub_field_total)
        annotation_info_dict = {}

        # Add annotations to the 
        #annotation_number = 0
        #if annotation_number < annotation_count:
        #    annotation_number =+ 1

        # Link each column in the mutation line with a field
        for field, mutation_info in zip(mut_fields, line):

            #mut_dict[field] = mutation_info
            
            # A list to capture all annotations and subfield data
            mut_annotation_list = []
            
            # For the field "INFO" store annotation information as a list of dicts
            if field == "INFO":

                # Annotation information is delimited by a "|"
                annotation_info = mutation_info.split("|")
                
                sub_field_counter = 0

                # Store the subfields as a list in each annotation object in the annotation dict
                for annotation_index in range(annotation_count):

                    #print("annotation_index " + str(annotation_index) + " out of " + str(annotation_count))

                    annotation_info_dict[annotation_index] = annotation_info[annotation_index*(sub_field_total):(annotation_index+1)*(sub_field_total-1)]
                
                #print(annotation_info_dict)
                #print("---------------------")

                
                # Add the annotation dictionary to the annotation list
                mut_annotation_list.append(annotation_info_dict)

            else:
                # Convert each mutation row into a dictionary with keys as fields mutation info as values
                mut_dict[field] = mutation_info
            
            # Add the annotation list to the mutation dictionary
            mut_dict["INFO"] = mut_annotation_list
        
        dict_list.append(mut_dict)
        #print(dict_list)

        # Only go through a limited number of rows for troubleshooting
        count += 1
        if count > stop_count:
            break

    # Export the mutation and their annotations as a json file
    output_json = str(output_folder) + "/" + str(output_prefix) + ".json"

    with open(output_json, 'w', encoding = 'utf-8') as output_json:
        json.dump(dict_list, output_json, indent=4)
    output_json.close()


    #print(dict_list)
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for civic vcf to csv convertor.')
    parser.add_argument('--input_vcf_file', '-i',
                        help='An absolute path to the vcf file from CIVIC')
    parser.add_argument('--schema', '-s',
                        help='A schema file of the subfields')                       
    parser.add_argument('--output_prefix', '-p',
                        help='A prefix for naming the output files')
    parser.add_argument('--output_folder', '-o',
                        help='An absolute path to the output folder')
    args = parser.parse_args()

    main(args.input_vcf_file, args.schema, args.output_prefix, args.output_folder)