'''
Input:
########
    * -i : A path to the ICGC .csv file
    * -m : A path to the folder containing mapping files
    * -d : The name of the doid mapping file
    * -e : The name of the ensp to uniprot accession mapping file
    * -o : A path to the output folder

Output:
########
    * A .csv file with mutation data formatted to the biomuta field structure

Usage:
########
    * python map_icgc.py -h

    *Gives a description of the neccessary commands

    * python map_icgc.py -i <path/input_file.vcf> -s <path/schema.json> -p output_prefix_name -o <path/>

    *Runs the script with the given csv file and outputs a csv file formatted for the final biomuta master file

'''

import argparse
from cmath import nan
import csv
import re
import pandas as pd



def main(icgc_csv, mapping_folder, doid_mapping_csv, enst_mapping_csv, output_folder):
    ##################################
    # Load the mapping files
    ##################################
    # Load in the TCGA mapping file to a mapping and cancer list
    doid_file_csv = mapping_folder + '/' + doid_mapping_csv
    enst_file_csv = mapping_folder + '/' + enst_mapping_csv
    
    with open(doid_file_csv, "r") as doid_mapping_handle:
        doid_mapping = csv.reader(doid_mapping_handle)
        # Skip the header
        next(doid_mapping)

        # Set up the mapping dictionary
        doid_mapping_dict = {}

        # Populate the mapping dict
        for row in doid_mapping:
            doid_mapping_dict[row[0]] = row[1]
        
        for key,value in doid_mapping_dict.items():
            re.sub(r'NA', '', value)
    
    # Set up a list of cancers to iterate through
    cancer_list = []
    for key, value in doid_mapping_dict.items():
        if value not in cancer_list:
            cancer_list.append(value)
    
    # Load the ENSP to uniprot mapping file.
    with open(enst_file_csv, "r") as enst_file_handle:
        enst_mapping = csv.reader(enst_file_handle, quoting=csv.QUOTE_ALL)
        # Skip the header.
        next(enst_mapping)

        # Set up the mapping file dictionary.
        ensp_mapping_dict = {}

        # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
        for row in enst_mapping:
            ensp_mapping_dict[row[2]] = row[0]
    
    columns_from_icgc = [
        'chr_id',
        'start_pos', 
        'ref_nt', 
        'alt_nt', 
        'transcript_name', 
        'Mutation genome position'
        ]
    
    print("Loading mutations from COSMIC tsv mutation file for fields: " + str(columns_from_cosmic))

    cosmic_df_iterator = pd.read_csv(cosmic_tsv, usecols=columns_from_cosmic, dtype=str, sep='\t',chunksize=1000000)

    for i, cosmic_df in enumerate(cosmic_df_iterator):

        print("Starting chunk " + str(i))

        # Map doid child to parent terms
        cosmic_df['do_name'] = cosmic_df['Primary site'].map(doid_mapping_dict)
    
        # Create new fields for reformatted data: ENST, genome location, AA mutation, nuceotide mutation
        cosmic_df['ENST'] = ''
        cosmic_df['chr_id'] = ''
        cosmic_df['start_pos'] = ''
        cosmic_df['end_pos'] = ''
        cosmic_df['ref_aa'] = ''
        cosmic_df['alt_aa'] = ''
        cosmic_df['aa_pos'] = ''
        cosmic_df['ref_nt'] = ''
        cosmic_df['alt_nt'] = ''
        cosmic_df.rename(columns = {'Sample name':'sample_name'},inplace=True)
    
        # Drop rows with missing information
        cosmic_df.dropna(subset=['Mutation AA','Mutation CDS','Mutation genome position'],inplace=True)
        total_rows = len(cosmic_df.index)
    
        # Create a new column with only ENST ID to be used for mapping. Also separate the AA notation, genome locations, and nucleotide change
        
        # Format nucleotide change and remove indels
        print('Formatting nucleotide information')
        cosmic_df['nucleotide_info'] = cosmic_df['Mutation CDS'].apply(nt_format)
        cosmic_df.dropna(subset=['nucleotide_info'],inplace=True)
        cosmic_df[['ref_nt','alt_nt']] = pd.DataFrame(cosmic_df['nucleotide_info'].tolist(), index=cosmic_df.index)
    
        # Format the amino acid change and position
        print('Formatting amino acid information')
        cosmic_df['amino_acid_info'] = cosmic_df['Mutation AA'].apply(aa_format)
        cosmic_df.dropna(subset=['amino_acid_info'],inplace=True)
        cosmic_df[['ref_aa','alt_aa','aa_pos']] = pd.DataFrame(cosmic_df['amino_acid_info'].tolist(), index=cosmic_df.index)
    
        # Format the genomic location
        print('Formatting genomic location information')
        cosmic_df['genome_location_info'] = cosmic_df['Mutation genome position'].apply(gen_location_format)
        cosmic_df.dropna(subset=['genome_location_info'],inplace=True)
        cosmic_df[['chr_id','start_pos','end_pos']] = pd.DataFrame(cosmic_df['genome_location_info'].tolist(), index=cosmic_df.index)
        
        # Map ENST symbol to uniprot accession
        print('Mapping ENST IDs to uniprot accession')
        cosmic_df['ENST'] = cosmic_df['Accession Number'].apply(lambda x : x.split('.')[0])
        cosmic_df['uniprotkb_canonical_ac'] = cosmic_df['ENST'].map(ensp_mapping_dict)
    
        cosmic_df.dropna(subset=['ref_nt', 'ref_aa', 'chr_id'],inplace=True)
     
        # Select and rename fields for integration with other sources
        final_fields = [
            'sample_name',
            'chr_id',
            'start_pos',
            'end_pos',
            'ref_nt',
            'alt_nt',
            'aa_pos',
            'ref_aa',
            'alt_aa',
            'do_name',
            'uniprotkb_canonical_ac'
        ]
    
        final_df = cosmic_df[final_fields]
    
        final_df['source'] = 'cosmic'
    
        # Remove rows that did not map to uniprot canonical transcripts and duplicates
        final_df.dropna(subset=['uniprotkb_canonical_ac'],inplace=True)
        final_df.drop_duplicates(keep='first',inplace=True)
        
        # How to handle the df chunks to process
        mode = 'w' if i == 0 else 'a'
        header = i == 0
        

        mapped_new_file_path = output_folder + "/cosmic_missense_biomuta_v5.csv"
        print("Adding processed data to " + mapped_new_file_path)
        final_df.to_csv(mapped_new_file_path, index = False, header=header, mode=mode)

        print("Chunk number " + str(i) + " completed")    
    
    
    
    
    ##################################
    # Load the civic csv file and map, then export
    ##################################
    icgc_df = pd.read_csv(icgc_csv, dtype=str)

    # Map doid child to parent terms
    civic_df['do_name'] = civic_df['CIViC Entity Disease'].map(doid_mapping_dict)
    civic_df['do_name'] = civic_df['do_name'].apply(convert_NA)
    
    # Create a column that removes the dot notation from the ENST IDs in civic data
    civic_df['sample_name'] = ''
    civic_df['ENST'] = ''
    civic_df['ref_aa'] = ''
    civic_df['alt_aa'] = ''
    civic_df['aa_pos'] = ''
    civic_df['source'] = 'civic'
    civic_df['end_pos'] = ''

    # Check for indels and remove
    civic_df['REF'] = civic_df['REF'].apply(remove_indels)

    # Format the amino acid change and position
    print('Formatting amino acid information')
    # amino acid changes to exclude
    civic_df['amino_acid_info'] = civic_df['CIViC Variant Name'].apply(aa_format)
    civic_df.dropna(subset=['amino_acid_info'],inplace=True)
    civic_df[['ref_aa','alt_aa','aa_pos']] = pd.DataFrame(civic_df['amino_acid_info'].tolist(), index=civic_df.index)
    
    # Create a new column with only ENST ID to be used for mapping and separate AA notation
    civic_df['ENST'] = civic_df['Feature'].apply(lambda x : x.split('.')[0])

    # Map ENST symbol to uniprot accession
    civic_df['uniprotkb_canonical_ac'] = civic_df['ENST'].map(ensp_mapping_dict)

    # Select and rename fields for integration with other sources
    final_fields = [
        'sample_name',
        '#CHROM',
        'POS',
        'end_pos',
        'REF',
        'ALT',
        'aa_pos',
        'ref_aa',
        'alt_aa',
        'do_name',
        'uniprotkb_canonical_ac',
        'source'
    ]

    final_df = civic_df[final_fields]

    # Name the fields in the exported file according tot he BioMuta convention
    final_df.rename(columns={
        '#CHROM': 'chr_id', 
        'POS': 'start_pos', 
        'REF': 'ref_nt', 
        'ALT': 'alt_nt'
    }, inplace=True)

    final_df['end_pos'] = final_df['start_pos']

    final_df.dropna(inplace=True)

    mapped_new_file_path = output_folder + "/civic_missense_biomuta_v5.csv"
    print("Exporting mapped file to " + mapped_new_file_path)
    final_df.to_csv(mapped_new_file_path, index = False)

###############################
# Functions for formatting data
###############################

# Format the amino acid infomation
def aa_format(aa_info):
    # Define exceptions and additonal information to remove
    aa_exceptions = ['FRAMESHIFT','RS','rs','MUTATION','fs','FS','DEL','c.','DUP','HOM']
    aa_clean_up = ['BCR-ABL_','PML-RARA_','EM4-ALK_','ALK_Fusion_','HIP1-ALK_','FIP1L1-PDGFRA_','CD74-ROS1_','ETV6-NTRK3_']
    exception_flag = 0
    for info in aa_clean_up:
        aa_info = re.sub(info, '', str(aa_info))   
    for exception in aa_exceptions:
        if re.search(exception,aa_info):
            exception_flag = 1
            aa_list = [nan,nan,nan]

    # Format the amino acid change
    if exception_flag == 0:
        aa_list = re.findall(r'[A-Z]',aa_info)
        # Account for nonsense mutations
        if re.search(r'\*',aa_info):
            stop_character = re.findall(r'\*',aa_info)
            aa_list.append(stop_character[0])
        
        aa_position = re.findall(r'\d+',aa_info)
        aa_list.append(aa_position[0])

    # Account for additional outlier cases
    if len(aa_list) != 3:
        return [nan,nan,nan]
    else:
        return aa_list
    
def remove_indels(nt_info):
    if len(nt_info) > 1:
        nt_info = nan
    
    return nt_info

def convert_NA(NA_value):
    if NA_value == 'NA':
        NA_value = nan
    
    return NA_value
    

                            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for civic mapping to doid and uniprot accessions.')
    parser.add_argument('--civic_csv', '-c',
                        help='An absolute path to the civic csv')
    parser.add_argument('--mapping_folder', '-m',
                        help='A path to the folder containing mapping files')                       
    parser.add_argument('--doid_mapping', '-d',
                        help='The name of the doid mapping file')
    parser.add_argument('--enst_mapping', '-e',
                        help='The name of the enst mapping file')
    parser.add_argument('--output_folder', '-o',
                        help='A path to the folder to export the mapped file')
    args = parser.parse_args()

    main(args.civic_csv, args.mapping_folder, args.doid_mapping, args.enst_mapping, args.output_folder)

#python map_civic_csv.py -c /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/civic/civic_mutations_Mar_01_2022.csv -m /mnt/c/Users/caule/github_general/biomuta/pipeline/convert_step2/mapping -d civic_doid_mapping.csv -e human_protein_transcriptlocus.csv -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled 