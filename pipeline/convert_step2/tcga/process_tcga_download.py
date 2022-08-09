'''
Input:
########
    * -i : A path to the ninput csv to reformat
    * -m : A path to the folder containing mapping files
    * -d : A path to the tcga study to doid mapping file
    * -e : A path to the ENSP to uniprot mapping file
    * -o : A path to the output folder


Output:
########
    * A data report comparing new AA sites to old AA sites for Biomuta

Usage:
########
    * python process_tcga_download.py -h

    *Gives a description of the neccessary commands

    * python process_tcga_download.py -i <path/input_file.vcf> -m <path/> -d <doid_mapping.csv -e <ensp_mapping.csv> -o <path/>

    *Runs the script with the given input tcga csv and output a formatted csv

'''

import argparse
import csv
import pandas as pd
import re

def main(input_csv, mapping_folder, doid_mapping_file, ensp_mapping_file, output_folder, ):
    ##################################
    # Load the mapping files
    ##################################
    # Load in the TCGA mapping file to a mapping and cancer list
    doid_file_csv = mapping_folder + '/' + doid_mapping_file
    ensp_file_csv = mapping_folder + '/' + ensp_mapping_file

    with open(doid_file_csv, "r") as mapping_handle:
        mapping_csv = csv.reader(mapping_handle, delimiter="\t")
        # Skip the header
        next(mapping_csv)

        # Set up the mapping dictionary
        mapping_dict = {}

        # Populate the mapping dict
        for row in mapping_csv:
            doid_term = str(row[0]) + " / " + str(row[1])
            mapping_dict[row[2]] = doid_term
    
    # Load the ENSP to uniprot mapping file.
    with open(ensp_file_csv, "r") as ensp_file_handle:
        ensp_file_csv = csv.reader(ensp_file_handle, quoting=csv.QUOTE_ALL)
        # Skip the header.
        next(ensp_file_csv)

        # Set up the mapping file dictionary.
        ensp_mapping_dict = {}

        # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
        for row in ensp_file_csv:
            ensp_mapping_dict[row[3]] = row[1]
    
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
    
    print("Loading mutations from TCGA csv mutation file for fields: " + str(field_list))

    icgc_df_iterator = pd.read_csv(input_csv, usecols=field_list, dtype=str, chunksize=1000000)

    for i, new_df in enumerate(icgc_df_iterator):

        print("Starting chunk " + str(i))

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
        
        print('Mapping and formatting complete. Adding to the final dataframe...')
        # Create the final dataframe to export
        final_fields = (
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
        )
        final_df = new_df.loc[:, final_fields]
        final_df.dropna(inplace=True)
        final_df.rename(columns={
            'case_barcode': 'sample_name',
            'Start_Position': 'start_pos',
            'End_Position': 'end_pos',
            'Reference_Allele': 'ref_nt',
            'Tumor_Seq_Allele1': 'alt_nt'
        }, inplace=True)
    
        final_df.drop_duplicates(keep='first',inplace=True)

        # How to handle the df chunks to process
        mode = 'w' if i == 0 else 'a'
        header = i == 0
    
        # Export the mapped new mutation data
        mapped_new_file_path = output_folder + "/tcga_missense_biomuta_v5.csv"
        final_df.to_csv(mapped_new_file_path, index = False, header=header, mode=mode)

        print("Chunk number " + str(i) + " completed") 

    print("Exporting processed data to " + mapped_new_file_path)

###############################
# Functions for formatting data
###############################

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
    parser = argparse.ArgumentParser(description='Commands for tcga reformatting and mapping to doid and uniprot accessions.')
    parser.add_argument('--input_csv', '-i',
                        help='An absolute path to the input tcga csv')
    parser.add_argument('--mapping_folder', '-m',
                        help='A path to the folder containing mapping files')                       
    parser.add_argument('--doid_mapping', '-d',
                        help='The name of the doid mapping file')
    parser.add_argument('--enst_mapping', '-e',
                        help='The name of the enst mapping file')
    parser.add_argument('--output_folder', '-o',
                        help='A path to the folder to export the mapped file')
    args = parser.parse_args()

    main(args.input_csv, args.mapping_folder, args.doid_mapping, args.enst_mapping, args.output_folder)


    # Example run:

#python process_tcga_download.py -i /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/tcga/TCGA_SNP_somatic_mutation_hg38.csv -m /mnt/c/Users/caule/github_general/biomuta/pipeline/convert_step2/mapping -d tcga_doid_mapping.csv -e human_protein_transcriptlocus.csv -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled 