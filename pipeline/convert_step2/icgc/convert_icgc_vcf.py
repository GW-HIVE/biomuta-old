'''
Input:
########
    * -i : A path to the ICGC .vcf file
    * -s : A schema file containing the field names in the annotations and to use for the output file
    * -o : A path to the output folder, where the transformed CSV data will go

Output:
########
    * A .csv file with mutation data where each row contains one mutation and one unique annotation

Usage:
########
    * python convert_icgc_vcf.py -h

    *Gives a description of the neccessary commands

    * python convert_icgc_vcf.py -i <path/input_file.vcf> -s <path/schema.json> -o <path/>

    *Runs the script with the given input vcf and schema json and outputs a csv file.

'''

import argparse
import csv
import re
import json

def main(input_vcf_file, schema, output_folder):
    '''
    Loads an input .vcf and converts rows to csv format 
    '''
    
    # A list to hold the vcf headers
    vcf_headers = []

    # A list to hold the mutation data
    data = []    

    # Separate the VCF headers and mutation data, store vcf headers for optional output
    print('Loading the input VCF...')
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
        
        print('Loaded ' + str(len(data)) + ' rows from input VCF')

    # Use the field schema file to create the output header
    print('Creating the output header...')
    with open(schema) as schema_file: 
        schema_template = json.load(schema_file)

    general_fields_schema = schema_template['schema'][0]['general_fields']
    consequence_subfield_schema = schema_template['schema'][0]['consequence_subfields']
    consequence_subfield_count = len(consequence_subfield_schema)
    occurrence_subfield_schema = schema_template['schema'][0]['occurrence_subfields']
    occurrence_subfield_count = len(occurrence_subfield_schema)
    additional_subfield_schema = schema_template['schema'][0]['additional_fields']

    out_header = []
    for field in general_fields_schema:
        out_header.append(field)
    for subfield in consequence_subfield_schema:
        out_header.append(subfield)
    for subfield in occurrence_subfield_schema:
        out_header.append(subfield)
    for field in additional_subfield_schema:
        out_header.append(field)

    output_csv = output_folder + '/' + 'icgc_converted_mutations.csv'
    
    # Used to show progress of script
    index = 0
    total_rows = len(data)
    
    # Write to the output file row by row
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


                        # Separate the list of all row information into sections
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
                        
                        # Process sections containing subfields that can have multiple annotations
                        consequence_dict = process_section(consequence_section_info, consequence_subfield_schema)
                        occurrence_dict = process_section(occurrence_section_info, occurrence_subfield_schema)
                        
                        # Add a row for each unique combination of a consequence and occurrence annotation combined with annotation-independent information
                        for consequence_number, consequence_info in consequence_dict.items():
    
                            for  occurrence_number, occurrence_info in occurrence_dict.items():

                                out_row = []
                                
                                # Add the row information so that it matches the order of the output header fields
                                for info in mut_general_info:
                                    out_row.append(info)
    
                                for subfield, value in consequence_info.items():
                                    out_row.append(value)
                                
                                for subfield, value in occurrence_info.items():
                                    out_row.append(value)
                                
                                for info in additional_info:
                                    out_row.append(info)

                                writer.writerow(out_row)
                             
# A function for processing each section of a VCF row
def process_section(section_info, section_subfields):
    
    # A dictionary to hold each unqique annotation given in a section
    annotation_info_dict = {}

    # Annotations are deliited by commas
    annotation_stop_index = str(section_info).count(',') + 1
    annotations = section_info.split(',')
    
    # Loop through a list of all subfields for all annotations
    # Separate annotations in the dictionary based on the number of the number of subfields per annotation for the section 
    for annotation_index in range(annotation_stop_index):
        if annotation_index == annotation_stop_index:
            break

        # Subfields in each annotation are delimited by pipes
        annotation_info = annotations[annotation_index].split('|')
        subfield_index = 0
        subfield_dict = {}

        for subfield in section_subfields:
            subfield_dict[subfield] = annotation_info[subfield_index]
            subfield_index +=1
        
        annotation_info_dict[annotation_index] = subfield_dict
    
    # Return a dictionary where the key is an annotation number and the values are key:value pairs of subfields and annotation speficic information
    return annotation_info_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for icgc vcf to csv convertor.')
    parser.add_argument('--input_vcf_file', '-i',
                        help='An absolute path to the vcf file from ICGC')
    parser.add_argument('--schema', '-s',
                        help='A schema file of the fields and section subfields')                       
    parser.add_argument('--output_folder', '-o',
                        help='An absolute path to the output folder')
    args = parser.parse_args()

    main(args.input_vcf_file, args.schema, args.output_folder)

    #python convert_icgc_vcf.py -i /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/icgc/icgc_missense_mutations_38.vcf -s /mnt/c/Users/caule/github_general/biomuta/pipeline/icgc/subfield_schema.json -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/icgc/