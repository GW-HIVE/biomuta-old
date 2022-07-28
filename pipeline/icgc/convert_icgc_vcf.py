'''
Input:
########
    * -i : A path to the ICGC .vcf file
    * -p : A prefix used for naming the output files
    * -o : A path to the output folder, where the data csv and vcf headers will go

Output:
########
    * A .csv file with mutation data and a .txt file with the .vcf headers

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
import re
import json
import numpy as np



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
                if re.search(r'##INFO=<ID=CONSEQUENCE', str(row)):
                    csq_sub_field_total = str(row).count('|') + 1
            
            # Pull out header for mutation data
            elif re.search(r'#CHROM', str(row)):
                
                mut_fields = row
            
            else:
                data.append(row)
    
    # A list of mutation dictionaries that will be converted into the output json
    dict_list = []

    # Load the annotation subfield schema and save the keys in a list
    with open(schema) as schema_file: 
        annotation_schema_template = json.load(schema_file)
    
    
    schema_key_list = []
    for key in annotation_schema_template:
        schema_key_list.append(key)

    # A list of mutation entries with additional commas in the info which breaks the annotation separator
    bad_row_list = []
    
    # Convert each line in the mutation data to json format
    for line in data: 

        # Check that any of the annotations for the mutation are non-silent
        if re.search(r'missense_variant', str(line)) or re.search(r'stop_gained', str(line)) or re.search(r'stop_lost', str(line)):

            # Separate the INFO column into 
    
    
            # Count the number of sub_fields in the info column
            # Formula compensates for the number of '|' being always 1 fewer for each annotation added
            sub_field_count = str(line).count('|') + int((str(line).count('|')/(csq_sub_field_total - 1 ))) 
            #print(sub_field_count)
    
            # A dictionary to hold all mutation data
            mut_dict = {}
    
            # A dictionary to hold the annotation subfields
            annotation_info_dict = {}
    
            # Link each column in the mutation line with a field
            for field, mutation_info in zip(mut_fields, line):
    
                #mut_dict[field] = mutation_info
                
                # A list to capture all annotations and subfield data
                mut_annotation_list = []
                mut_occurrence_list = []
                
                # For the field "INFO" store annotation information as a list of dicts
                if field == "INFO":
    
    
                    # Split the INFO field into the categories
                    info_cat_stop_index = str(mutation_info).count(';') + 1
    
                    for info_cat_index in range(info_cat_stop_index):
                        if info_cat_index == info_cat_stop_index:
                            break
                    
                    # Separate the info data by category
                    info_cat_all = mutation_info.count(';')
    
                    annotations_cat = info_cat_all[1]
    
                    # Annotation information is delimited by a "|"
                    #annotation_info = mutation_info.split('|')
                    annotation_stop_index = str(annotations_cat).count(',') + 1
                    annotations = annotations_cat.split(',')
    
                    # Check if row has additional commas that will break annotation separator
                    #sub_field_count = str(annotations_cat).count('|') + int((str(annotations_cat).count('|')/(csq_sub_field_total - 1 )))
                    
    
                    #if annotation_stop_index > (sub_field_count/csq_sub_field_total):
                    #    bad_row_list.append(line)
                    #    break
    
    
                    # Store the subfields as a list in each annotation object in the annotation dict
                    for annotation_index in range(annotation_stop_index):
    
                        if annotation_index == annotation_stop_index:
                            break
    
                        #print("Processing annotation " + str(annotation_index) + " out of " + str(range(annotation_stop_index)))
    
                        #annotation_info_dict[annotation_index] = annotations[annotation_index]
    
                        annotation_info = annotations[annotation_index].split('|')

                        # Check if the specific annotation is a non-silent mutation
                        if re.search(r'missense_variant', str(annotation_info)) or re.search(r'stop_gained', str(annotation_info)) or re.search(r'stop_lost', str(annotation_info)):

    
                            #annotation = annotation_info[annotation_index*(sub_field_total):(annotation_index+1)*(sub_field_total)]
        
        
                            #print(mutation_info)
                            #print(annotation_info)
                            #print(annotation)
        
                            sub_field = 0
        
                            sub_field_dict = {}
        
        
        
                            # Use the schema to assign each annotation information to a key
                            for key in annotation_schema_template:
                                sub_field_dict[key]= annotation_info[sub_field]
                                sub_field += 1
                                #print(key)
                                #print(sub_field)
                                #if sub_field == (sub_field_total): 
                                #    break
                                #else:
                                #    print(sub_field)
                                #    print(annotation[sub_field])
                                #    sub_field_dict[key]= annotation[sub_field]
                                
        
                            annotation_info_dict[annotation_index] = sub_field_dict
                            #annotation_info_dict[annotation_index] = annotation
                    
                    # Add the annotation dictionary to the annotation list
                    mut_annotation_list.append(annotation_info_dict)
    
                    occurrences_cat = info_cat_all[2]
    
                    occ_stop_index = str(occurrences_cat).count(',') + 1
                    occurrences = occurrences_cat.split(',')

                    occurrence_dict = {}

                    for occurrence_index in range(occ_stop_index):
                        
                        if occurrence_index == occ_stop_index:
                            break

                        occurence_info = occurrences.split('|')

                        for key in annotation_schema_template:
                            sub_field_dict[key]= occurence_info[sub_field]
                            sub_field += 1
                            occurrence_dict[occurrence_index] = sub_field_dict
                    
                    mut_occurrence_list.append(occurrence_dict)
                    

    
    
                    # Add the additional INFO categories as main fields
                    mut_dict["affected_donors"] = info_cat_all[3]
                    mut_dict["mutation"] = info_cat_all[4]
                    mut_dict["project_count"] = info_cat_all[5]
                    mut_dict["studies"] = info_cat_all[6]
                    mut_dict["tested_donors"] = info_cat_all[7]
    
                else:
                    # Convert each mutation row into a dictionary with keys as fields mutation info as values
                    mut_dict[field] = mutation_info
    
                
                # Add the annotation and occurrence list to the mutation dictionary
                mut_dict["INFO"] = mut_annotation_list
                mut_dict["Occurrence"] = mut_occurrence_list
            
            dict_list.append(mut_dict)
    
            # Only go through a limited number of rows for troubleshooting
            #count += 1
            #if count > stop_count:
            #    break

    
    # Create a a bad row error report
    if len(bad_row_list) > 0:
        print(str(len(bad_row_list)) + " entries contained additional commas could not be processed.")
        print("Entries are listed below:")
        print("----------------------------------")
        indexed_bad_row_list = ["%i: %s" % (index, value) for index, value in enumerate(bad_row_list)]
        formatted_bad_row_list = "\n\n".join(indexed_bad_row_list)
        print(formatted_bad_row_list)
        print("----------------------------------")
        print("Please only replace additional commas in these rows and rerun convertor. Leave commas that separate annotations unaltered.")
    
    # Export the mutation and their annotations as a json file
    output_json = str(output_folder) + "/" + str(output_prefix) + ".json"

    with open(output_json, 'w', encoding = 'utf-8') as output_json:
        json.dump(dict_list, output_json, indent=4)
    output_json.close()

    # Export a csv file with one annotation per row

    # Organize the csv headers
    csv_headers = []
    general_info_fields = []

    for field in mut_fields:
        csv_headers.append(field)
        general_info_fields.append(field)
    
    for key in schema_key_list:
        csv_headers.append(key)
    
    # The INFO field will be replace by the annotation information fields
    csv_headers.remove('INFO')
    general_info_fields.remove('INFO')

    # Create the csv file to populate with data
    output_csv = str(output_folder) + "/" + str(output_prefix) + ".csv"
    
    # Write the data dict to the csv
    with open(output_csv, 'w', encoding='utf-8') as output_csv: 
        writer = csv.writer(output_csv)
        writer.writerow(csv_headers)
        # Iterate through the mutation objects to fill the csv row

        #count = 0 
        #stop_count = 5
        
        # Iterate through each mutation
        for mutation in dict_list:
            #if count > stop_count:
            #    break
            #count += 1

            # Iterate through each field then find the INFO column
            for field, general_value in mutation.items():
                if field == "INFO":

                    # Iterate through each annotation
                    for annotation in mutation["INFO"]:

                        # For each annotation in the INFO field create a new row
                        for annotation_number, annotation_info in annotation.items():

                            # Create a new empty row
                            row = []

                            # Add general information
                            for field in general_info_fields:
                                row.append(mutation[field])
                                
                            # Add annotation specific information
                            for subfield, value in annotation_info.items():
                                row.append(value)
                    
                            # Check that each row has the same number of columns
                            if len(row) != len(csv_headers):
                                print("Bad row: does not match  number of fields")
                                print("Row length is " + str(len(row)) + " and should be " + str(len(csv_headers)))
                                print(row) 
                                continue
                            else:
                                # Write the row to the csv if it passes check
                                writer.writerow(row)
                            

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