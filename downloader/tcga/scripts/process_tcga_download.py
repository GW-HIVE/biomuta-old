'''
Input:
########
    * -n : A path to the newly updated Biomuta mutation list
    * -p : A path to the previous update of Biomuta's mutation list
    * -o : A path to the output folder, where the data report will be exported
    * -v : The version of the old biomuta to compare to the newly updated mutation list
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



def main(new_file, previous_version_file, tcga_mapping_file, uniprot_mapping_file, output_folder):
    '''
    '''
    ##################################
    # Load the mapping files
    ##################################
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
    
    # Set up a list of cancers to iterate through
    cancer_list = []
    for key, value in mapping_dict.items():
        if value not in cancer_list:
            cancer_list.append(value)
    
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
    
    ##################################
    # Process the new mutation file
    ##################################
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

    # Start a list of all genes present in the new dataset
    total_gene_list = []
    
    # A dictionary to hold every cancer type
    new_cancer_dictionary = {}

    # For the comparison, remove the rows not mapped to uniprot accession
    new_compare_df = new_df[new_df['uniprotkb_ac'].notna()]


    comparison_report_dict = {}

    total_mutation_counter = 0

    # The mutaiton comparison library is structured as:
    # Top level -- cancer type dict
    # Mid level -- list of gene dicts, one dict per gene in a cancer type
    # Bottom level -- list of mutations within each gene dict for that cancer type
    for cancer in cancer_list:

        #print("------------------------------------------")
        #print("Processing " + cancer)
        #print("------------------------------------------")

        
        gene_list = []
        gene_dict = {}
        comparison_details = {}

        # Start a counter for number of mutations per cancer
        mutation_counter = 0

        for index, row in new_compare_df.iterrows():

            # Check the cancer type for the row and if the gene symbol was mapped to a uniprot accession
            if row['doid_term'] == cancer:

                mutation_counter += 1

                # Add the gene to a gene list for that cancer
                if row['uniprotkb_ac'] not in gene_list:
                    gene_list.append(row['uniprotkb_ac'])
                    #print('Registering ' + str(row['uniprotkb_ac']) + ' for ' + cancer)

                    mut_list = [row['Amino_acids']]
                    
                    # Add the gene as a key in the gene dictionary for that cancer
                    gene_dict[row['uniprotkb_ac']] = mut_list

                    # Add the gene to the list for the whole dataset
                    if row['uniprotkb_ac'] not in total_gene_list:
                        total_gene_list.append(row['uniprotkb_ac'])
                    
                else:
                    # Add the mutation to the existing mutation list for that gene and cancer
                    gene_dict[row['uniprotkb_ac']].append(row['Amino_acids'])

                #print("Added " + str(row['Amino_acids']) + " to " + str(row['uniprotkb_ac']))

        # Add the mutation and gene lists to the all cancer dictionary
        new_cancer_dictionary[cancer] = gene_dict

        # Add the number of genes and mutations to the comparison report dictionary
        comparison_details['Mutations in new update'] = mutation_counter
        comparison_details['Genes in new update'] = len(gene_list)
        comparison_report_dict[cancer] = comparison_details

        # Add to the total mutation counter
        total_mutation_counter += mutation_counter
    
    # Add a section for all cancers combined to the comparison report
    comparison_details_total = {}
    comparison_details_total['Mutations in new update'] = total_mutation_counter
    comparison_details_total['Genes in new update'] = len(total_gene_list)
    comparison_report_dict['Totals'] = comparison_details_total

    ##################################
    # Process the old mutation file
    ##################################
    # Load the old file as a dataframe
    old_df = pd.read_csv(previous_version_file)

    # Start a list of all genes present in the new dataset
    total_old_gene_list = []
    
    # A dictionary to hold every cancer type
    old_cancer_dictionary = {}


    total_old_mutation_counter = 0

    # The mutaiton comparison library is structured as:
    # Top level -- cancer type dict
    # Mid level -- list of gene dicts, one dict per gene in a cancer type
    # Bottom level -- list of mutations within each gene dict for that cancer type
    for cancer in cancer_list:

        #print("------------------------------------------")
        #print("Processing " + cancer)
        #print("------------------------------------------")

        
        gene_list = []
        gene_dict = {}

        # Start a counter for number of mutations per cancer
        mutation_counter = 0

        for index, row in old_df.iterrows():

            # Check the cancer type for the row and if the gene symbol was mapped to a uniprot accession
            if row['do_name'] == cancer:

                mutation_counter += 1

                # Combine the two amino acid columns for ref and alt
                amino_acid_change = row['ref_aa'] + '/' + row['alt_aa']

                # Add the gene to a gene list for that cancer
                if row['xref_id'] not in gene_list:
                    gene_list.append(row['xref_id'])
                    #print('Registering ' + str(row['xref_id']) + ' for ' + cancer)

                    mut_list = [amino_acid_change]
                    
                    # Add the gene as a key in the gene dictionary for that cancer
                    gene_dict[row['xref_id']] = mut_list

                    # Add the gene to the list for the whole dataset
                    if row['xref_id'] not in total_old_gene_list:
                        total_old_gene_list.append(row['xref_id'])
                    
                else:
                    # Add the mutation to the existing mutation list for that gene and cancer
                    gene_dict[row['xref_id']].append(amino_acid_change)

                #print("Added " + str(amino_acid_change) + " to " + str(row['xref_id']))

        # Add the mutation and gene lists to the all cancer dictionary
        old_cancer_dictionary[cancer] = gene_dict

        # Add the number of genes and mutations to the comparison report dictionary
        #comparison_details["Mutations in old update"] = mutation_counter
        #comparison_details["Genes in old update"] = len(gene_list)
        comparison_report_dict[cancer]['Mutations in old update'] = mutation_counter
        comparison_report_dict[cancer]['Genes in old update'] = len(gene_list)

        # Add to the total mutation counter
        total_old_mutation_counter += mutation_counter
    
    # Add a section for all cancers combined to the comparison report
    comparison_details_total['Mutations in old update'] = total_old_mutation_counter
    comparison_details_total['Genes in old update'] = len(total_old_gene_list)
    comparison_report_dict['Totals'] = comparison_details_total



    
    ##################################
    # Compare new vs. old genes and mutation
    ##################################
    



    # Compare gene lists and add info to comparison report
    shared_genes = []
    missing_genes = []

    for gene in total_old_gene_list:
        if gene in total_gene_list:
            shared_genes.append(gene)
        else: 
            missing_genes.append(gene)

    comparison_report_dict['Totals']['Shared genes'] = len(shared_genes)
    comparison_report_dict['Totals']['Missing genes'] = len(missing_genes)

    # Compare mutations
    
    # Compare mutations per gene in each cancer
    for cancer in cancer_list:

        # Use the list of mutations from the old mutation list to find what mutations are missing in the new set
        for id, mutations in old_cancer_dictionary[cancer].items():
            for mutation in mutations:
                
                # Check if the gene is present in the new set
                if id in new_cancer_dictionary[cancer]:

                    if mutation in new_cancer_dictionary[cancer][id]:
                        print("Shared mutation")
                        print(id,mutation)

                    else: 
                        print("Missing mutation")
                        print(id,mutation)
                
                else:
                    #print(id + " not in new set")
                    pass



            


        
























    # Export the mutation lists and comparison report to json files
    output_new_json = str(output_folder) + "/new_TCGA_mutations_by_cancer.json"
    output_old_json = str(output_folder) + "/old_TCGA_mutations_by_cancer.json"
    comparison_report_json = str(output_folder) + "/TCGA_mutations_comparison_report.json"

    with open(output_new_json, 'w', encoding = 'utf-8') as output_new_json:
        json.dump(new_cancer_dictionary, output_new_json, indent=4)

    with open(output_old_json, 'w', encoding = 'utf-8') as output_old_json:
        json.dump(old_cancer_dictionary, output_old_json, indent=4)

    with open(comparison_report_json, 'w', encoding = 'utf-8') as comparison_report_json:
        json.dump(comparison_report_dict, comparison_report_json, indent=4)

                            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for civic vcf to csv convertor.')
    parser.add_argument('--new_file', '-n',
                        help='An absolute path to the new TCGA mutation download')
    parser.add_argument('--tcga_mapping_file', '-t',
                        help='An absolute path to the TCGA study to cancer mapping file')
    parser.add_argument('--previous_version_file', '-p',
                        help='An absolute path to the previous updates mutation file')
    parser.add_argument('--output_folder', '-o',
                        help='An absolute path to the output folder')
    parser.add_argument('--uniprot_mapping_file', '-u',
                        help='An absolute path to the uniprot accession mapping file')   
    args = parser.parse_args()

    main(args.new_file, args.previous_version_file, args.tcga_mapping_file, args.uniprot_mapping_file, args.output_folder)


    # Example run:

# python process_tcga_download.py -n /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/new_file_top.csv -p /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/old_mutations_top.csv -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/ -u mapping/uniprot_masterlist.csv -t mapping/TCGA_DOID_mapping_v4.0.csv