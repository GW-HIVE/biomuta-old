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
    


    # For the comparison, remove the rows not mapped to uniprot accession
    new_compare_df = new_df[new_df['uniprotkb_ac'].notna()]

    total_new_rows = len(new_compare_df.index)

    # The mutaiton comparison library is structured as:
    # Top level -- cancer type dict
    # Mid level -- list of gene dicts, one dict per gene in a cancer type
    # Bottom level -- list of mutations within each gene dict for that cancer type

    print("------------------------------------------")
    print("Processing new mutation file")
    print("------------------------------------------")
    

    # Need to set up blank dictionaries and lists here, then populate the master dictionary at each row instead of compiling the gene lists then adding them at the end

    # A dictionary to hold every cancer type
    new_cancer_dictionary = {}    
    comparison_report_dict = {}
    
    for index, row in new_compare_df.iterrows():

        print('Processing row ' + str(index) + ' out of ' + str(total_new_rows), end ='\r')
        #print('Processing row ' + str(index) + ' out of ' + str(total_new_rows))
        
        # Create the correct amino acid change formatting, skip outlier cases with different AA format listed (example 'X' instead of 'E/X')
        if re.search(r'/', str(row['Amino_acids'])):
            protein_position = str(row['Protein_position']).split('/')
            aa_position = protein_position[0]
            only_aa = str(row['Amino_acids']).split('/')
            amino_acid_change = only_aa[0] + aa_position + only_aa[1]
        else:
            continue
        

        # Start a counter for number of mutations per cancer
        for cancer in cancer_list:

            # Check the cancer type for the row and if the gene symbol was mapped to a uniprot accession
            if row['doid_term'] == cancer:

                # Create the cancer term in the gene dictionary
                if cancer in new_cancer_dictionary:
                    pass  
                else:
                    new_cancer_dictionary[cancer] = {}

                # Add the gene to a gene list for that cancer
                if row['uniprotkb_ac'] not in new_cancer_dictionary[cancer]:

                    # A list to hold the mutations for that gene
                    mut_list = [amino_acid_change]

                    # Assign the mutation list to the newly registered gene
                    new_cancer_dictionary[cancer][row['uniprotkb_ac']] = mut_list

                    # Add the gene to the list for the whole dataset
                    if row['uniprotkb_ac'] not in total_gene_list:
                        total_gene_list.append(row['uniprotkb_ac'])
                    
                else:
                    # Add the mutation to the existing mutation list for that gene and cancer
                    new_cancer_dictionary[cancer][row['uniprotkb_ac']].append(amino_acid_change)
                
                # Update or add cancer in the comparison dictionary
                if cancer in comparison_report_dict:
                    comparison_report_dict[cancer]['Mutations in new update'] += 1             

                else:
                    comparison_report_dict[cancer] = {}
                    comparison_report_dict[cancer]['Mutations in new update'] = 1
                    
                comparison_report_dict[cancer]['Genes in new update'] = len(new_cancer_dictionary[cancer].keys())
    
    # Add a section for all cancers combined to the comparison report
    comparison_details_total = {}
    comparison_details_total['Mutations in new update'] = total_new_rows
    comparison_details_total['Genes in new update'] = len(total_gene_list)
    comparison_report_dict['Totals'] = comparison_details_total


    output_new_json = str(output_folder) + "/new_TCGA_mutations_by_cancer.json"
    with open(output_new_json, 'w', encoding = 'utf-8') as output_new_json:
        json.dump(new_cancer_dictionary, output_new_json, indent=4)

    ##################################
    # Process the old mutation file
    ##################################
    # Load the old file as a dataframe
    old_df = pd.read_csv(previous_version_file)

    total_old_rows = len(old_df.index)

    # Start a list of all genes present in the new dataset
    total_old_gene_list = []

    # The mutaiton comparison library is structured as:
    # Top level -- cancer type dict
    # Mid level -- list of gene dicts, one dict per gene in a cancer type
    # Bottom level -- list of mutations within each gene dict for that cancer type
    
    print("------------------------------------------")
    print("Processing old mutation file")
    print("------------------------------------------")

    # A dictionary to hold every cancer type
    old_cancer_dictionary = {}

    absent_cancer_type_list = []

    for index, row in old_df.iterrows():

        print('Processing row ' + str(index) + ' out of ' + str(total_old_rows), end ='\r')

        for cancer in cancer_list:

            if row['do_name'] not in cancer_list:

                if row['do_name'] in absent_cancer_type_list:
                    pass
                else:
                    print('Mutations excluded for ' + str(row['do_name']) + ', this cancer type is not in the cancer type list.')
                    absent_cancer_type_list.append(row['do_name'])

            # Check the cancer type for the row and if the gene symbol was mapped to a uniprot accession
            if row['do_name'] == cancer:

               # Create the cancer term in the gene dictionary
                if cancer in old_cancer_dictionary:
                    pass  
                else:
                    old_cancer_dictionary[cancer] = {}

                # Combine the two amino acid columns for ref and alt
                amino_acid_change = str(row['ref_aa']) + str(row['aa_pos']) + str(row['alt_aa'])

                # Add the gene to a gene list for that cancer
                if row['xref_id'] not in old_cancer_dictionary[cancer]:

                    mut_list = [amino_acid_change]

                    # Assign the mutation list to the newly registered gene
                    old_cancer_dictionary[cancer][row['xref_id']] = mut_list

                    # Add the gene to the list for the whole dataset
                    if row['xref_id'] not in total_old_gene_list:
                        total_old_gene_list.append(row['xref_id'])
                    
                else:
                    # Add the mutation to the existing mutation list for that gene and cancer
                    old_cancer_dictionary[cancer][row['xref_id']].append(amino_acid_change)

                # Update or add cancer in the comparison dictionary
                if cancer in comparison_report_dict:
                    if 'Mutations in old update' in comparison_report_dict[cancer]:
                        comparison_report_dict[cancer]['Mutations in old update'] += 1
                    else:
                        comparison_report_dict[cancer]['Mutations in old update'] = 1

                else:
                    comparison_report_dict[cancer] = {}
                    comparison_report_dict[cancer]['Mutations in old update'] = 1
                    
                comparison_report_dict[cancer]['Genes in old update'] = len(old_cancer_dictionary[cancer].keys())
            
    
    # Add a section for all cancers combined to the comparison report
    comparison_details_total['Mutations in old update'] = total_old_rows
    comparison_details_total['Genes in old update'] = len(total_old_gene_list)
    comparison_report_dict['Totals'] = comparison_details_total



    
    ##################################
    # Compare new vs. old genes and mutation
    ##################################

    print("------------------------------------------")
    print("Comparing Mutation Sets")
    print("------------------------------------------")
    

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

    total_shared_mutation_counter = 0
    total_missing_mutations_counter = 0
    
    # Compare mutations per gene in each cancer
    for cancer in cancer_list:

        print("Comparing Mutations in " + cancer)

        shared_mutation_counter = 0
        missing_mutation_counter = 0

        # Check if cancer exists in the old mutation set
        if cancer in old_cancer_dictionary:

            # Use the list of mutations from the old mutation list to find what mutations are missing in the new set
            for id, mutations in old_cancer_dictionary[cancer].items():
    
                for mutation in mutations:
                    
                    # Check if cancer exists in new set
                    if cancer in new_cancer_dictionary:
    
                        # Check if the gene is present in the new set
                        if id in new_cancer_dictionary[cancer]:
        
                            if mutation in new_cancer_dictionary[cancer][id]:
                                shared_mutation_counter += 1
        
                            else: 
                                missing_mutation_counter += 1
                        
                        else:
                            missing_mutation_counter += 1
                    else:
                        missing_mutation_counter += 1
                    
            
            # Add the shared and missing mutation counts to the data report
            total_shared_mutation_counter += shared_mutation_counter
            comparison_report_dict[cancer]['Shared Mutations'] = shared_mutation_counter
            total_missing_mutations_counter += missing_mutation_counter
            comparison_report_dict[cancer]['Missing Mutations'] = missing_mutation_counter
    
    comparison_report_dict['Totals']['Shared Mutations'] = total_shared_mutation_counter
    comparison_report_dict['Totals']['Missing Mutations'] = total_missing_mutations_counter



    ##################################
    # Test genes
    ##################################

    print("------------------------------------------")
    print("Processing test genes")
    print("------------------------------------------")

    test_gene_list = ['P00533','Q13315','P01116','P42345','Q14118']

    comparison_report_dict['Test Genes'] = {}

    # Compare mutations between old and new sets for select genes
    for gene in test_gene_list:

        comparison_report_dict['Test Genes'][gene] = {}
        comparison_report_dict['Test Genes'][gene]['Shared Mutations'] = []
        comparison_report_dict['Test Genes'][gene]['Missing Mutations'] = []

        shared_gene_mutation_counter = 0
        missing_gene_mutation_counter = 0

        # Compare in each cancer
        for cancer in cancer_list:
            
            if cancer in old_cancer_dictionary:
    
                # Go through mutation list per gene
                for id, mutations in old_cancer_dictionary[cancer].items():
                    if id == gene:
                        for mutation in mutations:

                            # Do both cancer and gene exist in the new mutation set?
                            if cancer in new_cancer_dictionary:
                                if id in new_cancer_dictionary[cancer]:
    
                                    # After comparison add the mutation to eoither shared or missing 
                                    if mutation in new_cancer_dictionary[cancer][id]:
                                        shared_gene_mutation_counter += 1
                                        comparison_report_dict['Test Genes'][gene]['Shared Mutations'].append(mutation)
            
                                    else: 
                                        missing_gene_mutation_counter += 1
                                        comparison_report_dict['Test Genes'][gene]['Missing Mutations'].append(mutation)
                                
                                # Gene not found in new mutation set
                                else:
                                    missing_gene_mutation_counter += 1
                                    comparison_report_dict['Test Genes'][gene]['Missing Mutations'].append(mutation)

                            # Cancer not found in new mutation set
                            else:
                                missing_gene_mutation_counter += 1
                                comparison_report_dict['Test Genes'][gene]['Missing Mutations'].append(mutation)
    
            # Add the mutation counters to the comparison report
            comparison_report_dict['Test Genes'][gene]['Shared Mutation Count'] = shared_gene_mutation_counter
            comparison_report_dict['Test Genes'][gene]['Missing Mutation Count'] = missing_gene_mutation_counter

    print("------------------------------------------")
    print("Creating comparison report....")
    print("------------------------------------------")

    # Export the mutation lists and comparison report to json files
    #output_old_json = str(output_folder) + "/old_TCGA_mutations_by_cancer.json"
    comparison_report_json = str(output_folder) + "/TCGA_mutations_comparison_report.json"

    #with open(output_old_json, 'w', encoding = 'utf-8') as output_old_json:
    #    json.dump(old_cancer_dictionary, output_old_json, indent=4)

    print("------------------------------------------")
    print("All done! Comparison report writing to " + comparison_report_json)
    print("------------------------------------------")

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

#python process_tcga_download.py -n /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/tcga/TCGA_SNP_somatic_mutation_hg38.csv -p /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/tcga/human_protein_mutation_cancer.csv -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/ -u mapping/uniprot_masterlist.csv -t mapping/TCGA_DOID_mapping_v4.0.csv