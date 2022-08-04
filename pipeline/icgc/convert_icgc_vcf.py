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
            
            # Pull out header for mutation data
            elif re.search(r'#CHROM', str(row)):
                
                mut_fields = row
            
            else:
                data.append(row)

    # Create section schemas to be used for separating annotations
    with open(schema) as schema_file: 
        schema_template = json.load(schema_file)
    
    consequence_subfield_schema = schema_template['schema'][0]['consequence_subfields']
    consequence_subfield_count = len(consequence_subfield_schema)
    occurrence_subfield_schema = schema_template['schema'][0]['occurrence_subfields']
    occurrence_subfield_count = len(occurrence_subfield_schema)

    out_header = ['chr_id','id','start_pos','ref_nt','alt_nt','qual','filter']
    for subfield in consequence_subfield_schema:
        out_header.append(subfield)
    for subfield in occurrence_subfield_schema:
        out_header.append(subfield)

    # A list of mutation entries with additional commas in the info which breaks the annotation separator
    bad_row_list = []

    output_csv = output_folder + '/' + 'icgc_converted_mutations.csv'
    
    index = 0
    total_rows = len(data)

    with open(output_csv, 'w', encoding='utf-8') as output_csv: 
        writer = csv.writer(output_csv)
        writer.writerow(out_header)
    
        # Convert each line in the mutation data to json format
        for line in data: 
            
            index += 1

            print('Processing row ' + str(index) + ' out of ' + str(total_rows), end ='\r')
    
            # Check that any of the annotations for the mutation are non-silent
            if re.search(r'missense_variant', str(line)) or re.search(r'stop_gained', str(line)) or re.search(r'stop_lost', str(line)):
    
                # Separate the INFO column into 

                mut_general_info = []
        
                # Link each column in the mutation line with a field
                for field, mutation_info in zip(mut_fields, line):

                    if field != "INFO":
                        mut_general_info.append(mutation_info)
                    
                    # For the field "INFO" store annotation information as a list of dicts
                    if field == "INFO":
                        
                        # Split the INFO field into the categories
                        section_list = mutation_info.split(';')
    
                        if len(section_list) != 7:
                            print(line + ' contains irregular section count')
                            continue
    
                        consequence_section_info = section_list[0]

                        occurrence_section_info = section_list[1]

                        additional_info = []
                        affected_donors = section_list[2]
                        additional_info.append(affected_donors)
                        mutation = section_list[3]
                        additional_info.append(mutation)
                        project_info = section_list[4]
                        additional_info.append(project_info)
                        study_info = section_list[5]
                        additional_info.append(study_info)
                        tested_donors = section_list[6]
                        additional_info.append(tested_donors)
                            
                        consequence_dict = process_section(consequence_section_info, consequence_subfield_schema)
                        occurrence_dict = process_section(occurrence_section_info, occurrence_subfield_schema)
    
                        for consequence_number, consequence_info in consequence_dict.items():
    
                            for  occurrence_number, occurrence_info in occurrence_dict.items():

                                out_row = []

                                for info in mut_general_info:
                                    out_row.append(info)
    
                                for subfield, value in consequence_info.items():
                                    out_row.append(value)
                                
                                for subfield, value in occurrence_info.items():
                                    out_row.append(value)
                                
                                for info in additional_info:
                                    out_row.append(info)

                                writer.writerow(out_row)
    
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
    
    '''
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
    '''
                            

def process_section(section_info, section_subfields):

    subfields_per_annotation = len(section_subfields)

    subfield_total = str(section_info).count('|') + int((str(section_info).count('|')/(subfields_per_annotation - 1 )))

    annotation_info_dict = {}

    annotation_stop_index = str(section_info).count(',') + 1
    annotations = section_info.split(',')

    for annotation_index in range(annotation_stop_index):
        if annotation_index == annotation_stop_index:
            break
        annotation_info = annotations[annotation_index].split('|')
        subfield_index = 0
        subfield_dict = {}

        for subfield in section_subfields:
            subfield_dict[subfield] = annotation_info[subfield_index]
            subfield_index +=1
        
        annotation_info_dict[annotation_index] = subfield_dict
    
    return annotation_info_dict

    # Check if row has additional commas that will break annotation separator
    #if annotation_stop_index > (sub_field_count/sub_field_total):
        #bad_row_list.append(line)

    # Store the subfields as a list in each annotation object in the annotation dict
    #for annotation_index in range(annotation_stop_index
        #if annotation_index == annotation_stop_index:
        #    break
        #annotation_info = annotations[annotation_index].split('|')
        #sub_field = 0
        #sub_field_dict = 
        # Use the schema to assign each annotation information to a key
        #for key in annotation_schema_template:
        #    sub_field_dict[key]= annotation_info[sub_field]
        #    sub_field +=
        #annotation_info_dict[annotation_index] = sub_field_dict
    
    # Add the annotation dictionary to the annotation list
    #mut_annotation_list.append(annotation_info_dict)



    




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

    #python convert_icgc_vcf.py -i /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/icgc/EGFR_missense_icgc38.vcf -s /mnt/c/Users/caule/github_general/biomuta/pipeline/icgc/subfield_schema.json -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/icgc/