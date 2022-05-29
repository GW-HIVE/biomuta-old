'''
Input:
########
    * -n : A path to the newly updated Biomuta mutation list
    * -p : A path to the previous update of Biomuta's mutation list
    * -o : A path to the output folder, where the data report will be exported
    * -u : A path to the uniprot mapping file
    * -t : A path to the tcga study mapping file

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
import pandas as pd
import re
import json
import numpy as np



def main(new_file, previous_file, tcga_mapping_file, uniprot_mapping_file, output_folder):
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
            mapping_dict[row[2]] = doid_term
    
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
    new_df = pd.read_csv(new_file)

    # Map TCGA study names to doid terms
    new_df['doid_term'] = new_df['project_short_name'].map(mapping_dict)

    # Map HUGO symbol to uniprot id
    new_df['uniprotkb_ac'] = new_df['Hugo_Symbol'].map(uniprot_mapping_dict)


    # Export the mapped new mutation data
    mapped_new_file_path = output_folder + "/mapped_tcga_mutations.csv"
    print("Exporting mapped file to " + mapped_new_file_path)
    new_df.to_csv(mapped_new_file_path, index = False)

    # Set up a list of cancers to iterate through
    cancer_list = []
    for key, value in mapping_dict.items():
        if value not in cancer_list:
            cancer_list.append(value)

    
    # Start a list of all genes present in the dataset
    total_gene_list = []
    
    # A dictionary to hold every cancer type
    top_level_cancer_dict = {}

    # For the comparison, remove the rows not mapped to uniprot accession
    new_compare_df = new_df[new_df['uniprotkb_ac'].notna()]


    comparison_report_dict = {}

    total_mutation_counter = 0

    # The mutaiton comparison library is structured as:
    # Top level -- cancer type dict
    # Mid level -- list of gene dicts, one dict per gene in a cancer type
    # Bottom level -- list of mutations within each gene dict for that cancer type
    for cancer in cancer_list:

        print("------------------------------------------")
        print("Processing " + cancer)
        print("------------------------------------------")

        
        gene_list = []
        gene_dict = {}
        comparison_details_new = {}

        # Start a counter for number of mutations per cancer
        mutation_counter = 0

        for index, row in new_compare_df.iterrows():

            # Check the cancer type for the row and if the gene symbol was mapped to a uniprot accession
            if row['doid_term'] == cancer:

                mutation_counter += 1

                # Add the gene to a gene list for that cancer
                if row['uniprotkb_ac'] not in gene_list:
                    gene_list.append(row['uniprotkb_ac'])
                    print('Registering ' + str(row['uniprotkb_ac']) + ' for ' + cancer)

                    mut_list = [row['Amino_acids']]
                    
                    # Add the gene as a key in the gene dictionary for that cancer
                    gene_dict[row['uniprotkb_ac']] = mut_list

                    # Add the gene to the list for the whole dataset
                    if row['uniprotkb_ac'] not in total_gene_list:
                        total_gene_list.append(row['uniprotkb_ac'])
                    
                else:
                    # Add the mutation to the existing mutation list for that gene and cancer
                    gene_dict[row['uniprotkb_ac']].append(row['Amino_acids'])

                print("Added " + str(row['Amino_acids']) + " to " + str(row['uniprotkb_ac']))

        # Add the mutation and gene lists to the all cancer dictionary
        top_level_cancer_dict[cancer] = gene_dict

        # Add the number of genes and mutations to the comparison report dictionary
        comparison_details_new["Mutations in new update"] = mutation_counter
        comparison_details_new["Genes in new update"] = len(gene_list)
        comparison_report_dict[cancer] = comparison_details_new

        # Add to the total mutation counter
        total_mutation_counter += mutation_counter
    
    # Add a section for all cancers combined to the comparison report
    comparison_details_total = {}
    comparison_details_total["Mutations in new update"] = total_mutation_counter
    comparison_details_total["Genes in new update"] = len(total_gene_list)
    comparison_report_dict["Totals"] = comparison_details_total
    


    # Export the mutation lists and comparison report to json files
    output_json = str(output_folder) + "/TCGA_mutations_by_cancer.json"
    comparison_report_json = str(output_folder) + "/TCGA_mutations_comparison_report.json"

    with open(output_json, 'w', encoding = 'utf-8') as output_json:
        json.dump(top_level_cancer_dict, output_json, indent=4)
    output_json.close()

    with open(comparison_report_json, 'w', encoding = 'utf-8') as comparison_report_json:
        json.dump(comparison_report_dict, comparison_report_json, indent=4)
    comparison_report_json.close()


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
    #for cancer in cancer_list:

    #    print("Processing " + cancer)

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
    parser.add_argument('--new_file', '-n',
                        help='An absolute path to the new TCGA mutation download')
    parser.add_argument('--tcga_mapping_file', '-t',
                        help='An absolute path to the TCGA study to cancer mapping file')
    parser.add_argument('--previous_file', '-p',
                        help='An absolute path to the previous updates mutation file')
    parser.add_argument('--output_folder', '-o',
                        help='An absolute path to the output folder')
    parser.add_argument('--uniprot_mapping_file', '-u',
                        help='An absolute path to the uniprot accession mapping file')   
    args = parser.parse_args()

    main(args.new_file, args.previous_file, args.tcga_mapping_file, args.uniprot_mapping_file, args.output_folder)


    # Example run:

# python process_tcga_download.py -n /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/new_file_top.csv -p /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/old_file_example.csv -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/ -u mapping/uniprot_masterlist.csv -t mapping/TCGA_DOID_mapping_v4.0.csv