'''
Input:
########
    * -i : A path to the CIVIC .vcf file
    * -p : A prefix used for naming the output files
    * -o : A path to the output folder, where the mutation data csv will go
Output:
########
    * A .csv file with mutation data

Usage:
########
    * python convert_civic_vcf.py -h

    *Gives a description of the neccessary commands

    * python convert_civic_vcf.py -i <path/input_file.vcf> -s <path/schema.json> -o <path/>

    *Runs the script with the given input vcf and outputs a csv file.

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

    # Separate the VCF headers and mutation data
    print('Loading the input VCF...')
    with open(input_vcf_file, 'r') as vcf:
        
        reader = csv.reader(vcf, delimiter="\t")
        
        # Go through the VCF by row
        for row in reader:

            # Pull out VCF headers
            if re.search(r'##', str(row)):
                vcf_headers.append(row)
                
                # Record number of subfields in the INFO column
                if re.search(r'##INFO=<ID=CSQ', str(row)):
                    sub_field_total = str(row).count('|') + 1
            
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
    annotation_subfield_schema = schema_template['schema'][0]['annotation_subfields']

    out_header = []
    for field in general_fields_schema:
        out_header.append(field)
    for subfield in annotation_subfield_schema:
        out_header.append(subfield)

    output_csv = output_folder + '/' + 'civic_converted_mutations.csv'

    # Used to show progress of script
    index = 0
    total_rows = len(data)

    # A list of variant classifications to include in the final output
    variant_type_list = [
        'Missense_Variant',
        'Stop_Gained',
        'Stop_Lost',
        'Start_Gained',
        'Start_Lost'
    ]

    # Write to the output file row by row
    with open(output_csv, 'w', encoding='utf-8') as output_csv: 
        writer = csv.writer(output_csv)
        writer.writerow(out_header)
    
        # Convert each line in the mutation data to json format
        for line in data: 
            
            index += 1

            print('Processing row ' + str(index) + ' out of ' + str(total_rows), end ='\r')

            # Check that the mutation is a variant type to include int he final data
            variant_type_flag = 0
            for variant_type in variant_type_list:
                if re.search(variant_type, str(line)):
                    variant_type_flag = 1
    
            # Check that any of the annotations for the mutation are non-silent
            if variant_type_flag == 1:
    
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
    
                        if len(section_list) != 3:
                            print(line + ' contains irregular section count, should be 3 (sections separated by '';'')')
                            continue


                        # Separate the list of all row information into sections
                        gene_section_info = section_list[0]
                        variant_section_info = section_list[1]
                        annotation_section_info = section_list[2]
                        
                        # Process sections containing subfields that can have multiple annotations
                        gene_section_processed = process_gene_section(gene_section_info)
                        variant_section_processed = process_variant_section(variant_section_info) 
                        annotation_dict = process_annotation_section(annotation_section_info, annotation_subfield_schema)
                        if annotation_dict == None:
                            print("Skipping line for: ")
                            print(mut_general_info)
                            continue


                        # Add a row for each unique combination of a consequence and occurrence annotation combined with annotation-independent information
                        for annotation_number, annotation_info in annotation_dict.items():

                            out_row = []
                            
                            # Add the row information so that it matches the order of the output header fields
                            for info in mut_general_info:
                                out_row.append(info)

                            out_row.append(gene_section_processed)
                            
                            out_row.append(variant_section_processed)
                            
                            for subfield, info in annotation_info.items():
                                out_row.append(info)
                            writer.writerow(out_row)

# A function to remove extra charcaters in the gene section
def process_gene_section(section_info):
    processed_info = re.sub(r'GN=', '', section_info)
    return processed_info

# A function to remove extra charcaters in the amino acid variant section
def process_variant_section(section_info):
    processed_info = re.sub(r'VT=', '', section_info)
    return processed_info

# A function for processing each section of a VCF row
def process_annotation_section(section_info, section_subfields):

    subfield_count = len(section_subfields)
    
    # A dictionary to hold each unqique annotation given in a section
    annotation_info_dict = {}

    # Annotations are delimited by commas
    annotation_stop_index = str(section_info).count(',') + 1
    annotations = section_info.split(',')
    
    # Loop through a list of all subfields for all annotations
    # Separate annotations in the dictionary based on the number of the number of subfields per annotation for the section 
    for annotation_index in range(annotation_stop_index):
        if annotation_index == annotation_stop_index:
            break

        # Subfields in each annotation are delimited by pipes
        annotation_info = annotations[annotation_index].split('|')
        annotation_info_count = len(annotation_info)

        if annotation_info_count != subfield_count:
            print("Subfield count in annotation does not match expected subfield count. Inspect the following info and remove additional commas (likely between DOID terms)")
            print(section_info)
            return

        subfield_index = 0
        subfield_dict = {}

        for subfield in section_subfields:
            subfield_dict[subfield] = annotation_info[subfield_index]
            subfield_index +=1
        
        annotation_info_dict[annotation_index] = subfield_dict
    
    # Return a dictionary where the key is an annotation number and the values are key:value pairs of subfields and annotation speficic information
    return annotation_info_dict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for civic vcf to csv convertor.')
    parser.add_argument('--input_vcf_file', '-i',
                        help='An absolute path to the vcf file from CIVIC')
    parser.add_argument('--schema', '-s',
                        help='A schema file of the subfields')                       
    parser.add_argument('--output_folder', '-o',
                        help='An absolute path to the output folder')
    args = parser.parse_args()

    main(args.input_vcf_file, args.schema, args.output_folder)


# Example Run:
# python convert_civic_vcf.py -i /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/civic/01-Aug-2022-civic_accepted_hg38.vcf -s subfield_schema.json -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/civic