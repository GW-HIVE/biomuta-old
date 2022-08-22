'''
Input:
########
    * -c : A path to the cosmic tsv mutation file
    * -m : A path to the folder containing mapping files
    * -d : The name of the doid to cosmic cancer type mapping file
    * -e : The name of the enst to uniprot accession mapping file
    * -o : A path to the the folder to export the final mapped mutations


Output:
########
    * A mutation file with COSMIC mutations mapped to doid terms and uniprot accessions

Usage:
########
    * map_cosmic_tsv -h

    *Gives a description of the neccessary commands

    * python map_cosmic_tsv.py -c <path/cosmic_file_name.tsv> -m <path/mapping_folder> -d <doid_mapping_file_name> -e <enst_mapping_file_name> -o <path/output_folder>

    *Runs the script with the given input tsv and outputs a csv with Biomuta formatting.

'''

import argparse
from cmath import nan
import csv
import pandas as pd
import re

def main(cosmic_tsv, mapping_folder, doid_mapping_csv, enst_mapping_csv, output_folder):
    ##################################
    # Load the mapping files
    ##################################
    doid_file_csv = mapping_folder + '/' + doid_mapping_csv
    enst_file_csv = mapping_folder + '/' + enst_mapping_csv

    print("Loading DOID mapping file...")
    with open(doid_file_csv, "r") as doid_mapping_handle:
        doid_mapping = csv.reader(doid_mapping_handle)
        # Skip the header
        next(doid_mapping)

        # Set up the mapping dictionary
        doid_mapping_dict = {}

        # Populate the mapping dict
        for row in doid_mapping:
            doid_mapping_dict[row[0]] = row[1]
    
    # Load the ENST to uniprot mapping file.
    print("Loading ENST mapping file...")
    with open(enst_file_csv, "r") as enst_file_handle:
        enst_mapping = csv.reader(enst_file_handle, quoting=csv.QUOTE_ALL)
        # Skip the header.
        next(enst_mapping)

        # Set up the mapping file dictionary.
        ensp_mapping_dict = {}

        # Populate the mapping dictionary with keys as ensg IDs and values as the gene symbol.
        for row in enst_mapping:
            ensp_mapping_dict[row[2]] = row[1]
    
    ##################################
    # Load the cosmic tsv file and map, then export
    ##################################
    columns_from_cosmic = [
        'Accession Number',
        'Sample name', 
        'Primary site', 
        'Mutation CDS', 
        'Mutation AA', 
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
        final_fields = (
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
        )
    
        final_df = cosmic_df.loc[:, final_fields]
    
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


###############################
# Functions for formatting data
###############################

# Format the nucleotide information
def nt_format(nt_info):
    # Separate the nucleotide change
    nuc_list = re.findall(r'[A-Z]',nt_info)
    # Remove indels
    if len(nuc_list) != 2:
        #print('Excluding entry: ' + str(nt_info))
        #return pd.Series([nan,nan])
        nuc_list = [nan,nan]
    else:
        #return pd.Series(nuc_list, index=['ref_nt', 'alt_nt'])
        return (nuc_list)

# Format the amino acid infomation
def aa_format(aa_info):
    aa_list = re.findall(r'[A-Z\*]',aa_info)
    aa_position = re.findall(r'\d+',aa_info)
    aa_list.append(aa_position[0])
    if len(aa_list) != 3:
        return [nan,nan,nan]
    else:
        return aa_list

# Format the genomic location
def gen_location_format(gl_info):
    genome_info = str(gl_info).split(':')
    chr_id = genome_info[0]
    gen_loc_list = [chr_id]

    location_range = genome_info[1]
    location_positions = str(location_range).split('-')
    start_pos = location_positions[0]
    end_pos = location_positions[1]
    gen_loc_list.append(start_pos)
    gen_loc_list.append(end_pos)
    if len(gen_loc_list) != 3:
        return [nan,nan,nan] 
    else:
        return gen_loc_list


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for civic mapping to doid and uniprot accessions.')
    parser.add_argument('--cosmic_tsv', '-c',
                        help='An absolute path to the cosmic tsv')
    parser.add_argument('--mapping_folder', '-m',
                        help='A path to the folder containing mapping files')                       
    parser.add_argument('--doid_mapping', '-d',
                        help='The name of the doid mapping file')
    parser.add_argument('--enst_mapping', '-e',
                        help='The name of the enst mapping file')
    parser.add_argument('--output_folder', '-o',
                        help='A path to the folder to export the mapped file')
    args = parser.parse_args()

    main(args.cosmic_tsv, args.mapping_folder, args.doid_mapping, args.enst_mapping, args.output_folder)
        
#python map_cosmic_tsv.py -c /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/downloads/cosmic/Cosmic_SNPs_June_2022.tsv -m /mnt/c/Users/caule/github_general/biomuta/pipeline/convert_step2/mapping -d cosmic_doid_mapping.csv -e human_protein_transcriptlocus.csv -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled