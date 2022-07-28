'''
Input:
########
    * -n : A path to the newly updated Biomuta mutation list
    * -c : Should the new mutation sbe compared to an older mutation file?
    * -p : A path to the previous update of Biomuta's mutation list
    * -o : A path to the output folder, where the data report will be exported
    * -v : The version of the old biomuta to compare to the newly updated mutation list
    * -e : A path to the ENSP to uniprot mapping file
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
import pandas as pd
import re
import json
import numpy as np



def main(new_file, comparison, icgc_mutations, previous_version_file, tcga_mapping_file, ensp_mapping_file, output_folder, ):
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
    
    # Load the ENSP to uniprot mapping file.
    with open(ensp_mapping_file, "r") as ensp_file_handle:
        ensp_file_csv = csv.reader(ensp_file_handle, quoting=csv.QUOTE_ALL)
        # Skip the header.
        next(ensp_file_csv)

        # Set up the mapping file dictionary.
        ensp_mapping_dict = {}

        # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
        for row in ensp_file_csv:
            ensp_mapping_dict[row[3]] = row[0]
    
    ##################################
    # Process the new mutation file
    ##################################
    # Load the new file as a dataframe, process all values as strings to speed up conversion to dataframe
    field_list = [
        'project_short_name',
        'case_barcode',
        'Chromosome',
        'Start_Position',
        'End_Position',
        'Reference_Allele',
        'Tumor_Seq_Allele1',
        'Amino_acids',
        'Protein_position',
        'ENSP'
    ]
    
    print('Loading the mutations as a dataframe...')
    new_df = pd.read_csv(new_file, dtype=str, usecols=field_list)
    new_df.dropna(inplace=True)

    print('Mapping doid terms and uniprot accessions...')
    # Map TCGA study names to doid terms
    new_df['do_name'] = new_df['project_short_name'].map(mapping_dict)

    # Map HUGO symbol to uniprot id
    new_df['uniprotkb_canonical_ac'] = new_df['ENSP'].map(ensp_mapping_dict)

    # Format the amino acid ref, alt, and position columns, also the chromosome id
    print('Formatting amino acid info...')
    new_df['aa_pos'] = new_df['Protein_position'].apply(aa_format_position)
    new_df['aa_info'] = new_df['Amino_acids'].apply(aa_format_info)
    new_df['ref_aa'] = ''
    new_df['alt_aa'] = ''
    new_df[['ref_aa','alt_aa']] = pd.DataFrame(new_df['aa_info'].tolist(), index=new_df.index)
    new_df['chr_id'] = new_df['Chromosome'].apply(remove_chr)

    # Create a column for the source
    new_df['source'] = 'tcga'
    
    print('Mapping and formatting complete. Creating the final dataframe...')
    # Create the final dataframe to export
    final_fields = [
        'case_barcode',
        'chr_id',
        'Start_Position',
        'End_Position',
        'Reference_Allele',
        'Tumor_Seq_Allele1',
        'aa_pos',
        'ref_aa',
        'alt_aa',
        'do_name',
        'uniprotkb_canonical_ac',
        'source'
    ]
    final_df = new_df[final_fields]
    final_df.dropna(inplace=True)
    final_df.rename(columns={
        'case_barcode': 'sample_name',
        'Start_Position': 'start_pos',
        'End_Position': 'end_pos',
        'Reference_Allele': 'ref_nt',
        'Tumor_Seq_Allele1': 'alt_nt'
    }, inplace=True)

    final_df.drop_duplicates(keep='first',inplace=True)

    # Export the mapped new mutation data
    mapped_new_file_path = output_folder + "/tcga_missense_biomuta_v5.csv"
    print("Exporting mapped file to " + mapped_new_file_path)
    final_df.to_csv(mapped_new_file_path, index = False)

    # Start a list of all genes present in the new dataset
    total_gene_list = []

    # Check if mutation comparison should be performed
    if comparison == True:
        pass
    else:
        print('No comparison performed, all done!')
        exit()

    # For the comparison, remove the rows not mapped to uniprot accession
    new_compare_df = new_df[new_df['uniprotkb_canonical_ac'].notna()]

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
    

    ##################################
    # Add the ICGC mutations
    ##################################
    # Load the ICGC mutations and add to TCGA mutations if it was specified
    if icgc_mutations != 'no_icgc':

        print('Processing ICGC mutations...')
        
        # Process all the columns as strings to speed up dataset conversion to a dataframe
        icgc_df = pd.read_csv(icgc_mutations, dtype=str)
        
        total_icgc_rows = len(icgc_df.index)

        for index, row in icgc_df.iterrows():

            print('Processing row ' + str(index) + ' out of ' + str(total_icgc_rows), end ='\r')

            #Create the Amino Acid column 
            amino_acid_change = row['ref_aa'] + str(row['uniprot_pos']) + row['alt_aa']

            #Take the uniprot protein id without the isoform
            uniprot_info = str(row['uniprot_ac']).split('-')
            uniprot_id = uniprot_info[0]

            for cancer in cancer_list:

                if row['do_name'] == cancer:

                    # Create the cancer term in the gene dictionary
                    if cancer in new_cancer_dictionary:
                        pass  
                    else:
                        new_cancer_dictionary[cancer] = {}
    
                # Add the gene to a gene list for that cancer
                if uniprot_id not in new_cancer_dictionary[cancer]:

                    # A list to hold the mutations for that gene
                    mut_list = [amino_acid_change]

                    # Assign the mutation list to the newly registered gene
                    new_cancer_dictionary[cancer][uniprot_id] = mut_list

                    # Add the gene to the list for the whole dataset
                    if uniprot_id not in total_gene_list:
                        total_gene_list.append(uniprot_id)
                    
                else:
                    # Add the mutation to the existing mutation list for that gene and cancer
                    new_cancer_dictionary[cancer][uniprot_id].append(amino_acid_change)
                
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

        # For processing only TCGA data
        #if row['data_source'] == 'tcga':
        #    pass
        #else:
        #    continue

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

        shared_gene_mutation_counter = 0
        missing_gene_mutation_counter = 0

        # Compare in each cancer
        for cancer in cancer_list:



            
            if cancer in old_cancer_dictionary:
    
                # Go through mutation list per gene
                for id, mutations in old_cancer_dictionary[cancer].items():
                    if id == gene:
                        for mutation in mutations:
                            
                            # Set up lists per gene and per cancer in the test genes section
                            if gene not in comparison_report_dict['Test Genes']:
                                comparison_report_dict['Test Genes'][gene] = {}
                                
                            if cancer not in comparison_report_dict['Test Genes'][gene]:
                                comparison_report_dict['Test Genes'][gene][cancer] = {}
                                comparison_report_dict['Test Genes'][gene][cancer]['Shared Mutations'] = []
                                comparison_report_dict['Test Genes'][gene][cancer]['Missing Mutations'] = []

                            # Do both cancer and gene exist in the new mutation set?
                            if cancer in new_cancer_dictionary:
                                if id in new_cancer_dictionary[cancer]:
    
                                    # After comparison add the mutation to eoither shared or missing 
                                    if mutation in new_cancer_dictionary[cancer][id]:
                                        shared_gene_mutation_counter += 1
                                        comparison_report_dict['Test Genes'][gene][cancer]['Shared Mutations'].append(mutation)
            
                                    else: 
                                        missing_gene_mutation_counter += 1
                                        comparison_report_dict['Test Genes'][gene][cancer]['Missing Mutations'].append(mutation)
                                
                                # Gene not found in new mutation set
                                else:
                                    missing_gene_mutation_counter += 1
                                    comparison_report_dict['Test Genes'][gene][cancer]['Missing Mutations'].append(mutation)

                            # Cancer not found in new mutation set
                            else:
                                missing_gene_mutation_counter += 1
                                comparison_report_dict['Test Genes'][gene][cancer]['Missing Mutations'].append(mutation)
    
                        # Add the mutation counters to the comparison report
                        comparison_report_dict['Test Genes'][gene][cancer]['Shared Mutation Count'] = shared_gene_mutation_counter
                        comparison_report_dict['Test Genes'][gene][cancer]['Missing Mutation Count'] = missing_gene_mutation_counter

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

# Function to format the chromosome id
def remove_chr(chr_info):
    chr_info_update = re.sub(r'chr', '', chr_info)
    return chr_info_update

# Functions to format the amino acid change
def aa_format_position(aa_pos): 
    aa_pos_list = aa_pos.split('/')
    return aa_pos_list[0]

def aa_format_info(aa_info):
    aa_info_list = aa_info.split('/')
    return aa_info_list


                            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for civic vcf to csv convertor.')
    parser.add_argument('--new_file', '-n',
                        help='An absolute path to the new TCGA mutation download')
    parser.add_argument('--comparison', '-c',
                        help='Either yes or no to perform the mutation comparison', nargs='?', default=False)
    parser.add_argument('--icgc_mutations', '-i',
                        help='An absolute path to the ICGC mutation file', nargs='?', default='no_icgc')
    parser.add_argument('--tcga_mapping_file', '-t',
                        help='An absolute path to the TCGA study to cancer mapping file')
    parser.add_argument('--previous_version_file', '-p',
                        help='An absolute path to the previous updates mutation file', nargs='?', default='no_previous')
    parser.add_argument('--output_folder', '-o',
                        help='An absolute path to the output folder')
    parser.add_argument('--ensp_mapping_file', '-e',
                        help='An absolute path to the ensp to ensp to uniprot accession mapping file')   
    args = parser.parse_args()

    main(args.new_file, args.comparison, args.icgc_mutations, args.previous_version_file, args.tcga_mapping_file, args.ensp_mapping_file, args.output_folder)


    # Example run:

#python process_tcga_download.py -n /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/tcga/TCGA_SNP_somatic_mutation_hg38.csv -i /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/icgc/biomuta_3.0/biomuta_icgc.csv -p /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/tcga/human_protein_mutation_cancer.csv -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/ -u mapping/uniprot_masterlist.csv -t mapping/TCGA_DOID_mapping_v4.0.csv

#python process_tcga_download.py -n /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/tcga/TCGA_SNP_somatic_mutation_hg38.csv -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/ -e /mnt/c/Users/caule/github_general/biomuta/pipeline/convert_step2/mapping/human_protein_transcriptlocus.csv -t /mnt/c/Users/caule/github_general/biomuta/pipeline/convert_step2/mapping/TCGA_DOID_mapping_v4.0.csv