'''
Input:
########
    * -b : The final biomuta mutation file to extract mutations from
    * -a : The name of the uniprot accession for mutation extraction
    * -g : The path to the folder containing glycsy
    * -o : The output folder to export mutations for selected accessions

Output:
########
    * A csv file per gene, containing mutations extracted from the main biomuta mutation file

Usage:
########
    * python combine_csv.py -h

    *Gives a description of the neccessary commands

    * python combine_csv.py -b <path/final_biomuta_file.csv> -a <some_uniprot_accession> -o <path/>

    *Runs the script with the final biomuta mutation file and extracts mutations per accession given

'''

import argparse
from locale import D_T_FMT
import pandas as pd
import glob
import re



def main(biomuta_csv, accession, glyco_file_folder, output_folder):

    print('Using accession ' + accession + ' for QC')

    # Set up a dataframe to hold on mutations for the accession
    accession_df = pd.DataFrame()

    # Find all glycosylation datasets
    glyco_file_list = []
    glyco_file_path = glyco_file_folder + '*.csv'
    for file in glob.glob(glyco_file_path):
        glyco_file_list.append(file)

    
    # Set up dataframes to hold the glycosylation comparison
    for file in glyco_file_list:
        if re.search(r'loss', file):
            print('Glycoloss file: ' + file)
            glyco_loss_df = pd.read_csv(file, dtype = str)
            print('Rows in file: ' + str(len(glyco_loss_df.index)))
        elif re.search(r'effect', file):
            print('Glycoeffect file ' + file)
            glyco_effect_df = pd.read_csv(file, dtype=str)
            print('Rows in file: ' + str(len(glyco_effect_df.index)))
        else:
            print('Glycosylation dataset ' + file + ' not used for comparison')
    
    glyco_loss_overlap_df = pd.DataFrame()
    glyco_effect_overlap_df = pd.DataFrame()
    glyco_loss_unique_df = pd.DataFrame()
    glyco_effect_unique_df = pd.DataFrame()
    acc_glyco_loss_overlap_df = pd.DataFrame()
    acc_glyco_effect_overlap_df = pd.DataFrame()

    
    # Set up the chunk size for better memory usage and script speed
    biomuta_df_iterator = pd.read_csv(biomuta_csv, dtype=str, chunksize=1000000)
    
    # Add mutations to the dataframes using chunks of the master biomuta file
    for i, biomuta_df_chunk in enumerate(biomuta_df_iterator):

        print("Processing chunk number " + str(i+1))
        
        # Only gather mutations for the chosen accession
        acc_chunk_df = biomuta_df_chunk[biomuta_df_chunk['uniprotkb_canonical_ac'] == accession]
        
        # Run the glycosylation dataset comparison and add mutations to glycoloss or glycoeffect for the accession
        glyco_loss_overlap_df = pd.concat([glyco_loss_overlap_df, compare_shared_glyco_file(biomuta_df_chunk, glyco_loss_df)])
        glyco_effect_overlap_df = pd.concat([glyco_effect_overlap_df, compare_shared_glyco_file(biomuta_df_chunk, glyco_effect_df)])
        acc_glyco_loss_overlap_df = pd.concat([glyco_loss_overlap_df, compare_shared_glyco_file(acc_chunk_df, glyco_loss_df)])
        acc_glyco_effect_overlap_df = pd.concat([glyco_effect_overlap_df, compare_shared_glyco_file(acc_chunk_df, glyco_effect_df)])        

        # Glyco mutations not covered by the biomuta mutations
        glyco_loss_unique_df = pd.concat([glyco_loss_unique_df, compare_unique_glyco_file(biomuta_df_chunk, glyco_loss_df)])
        glyco_effect_unique_df = pd.concat([glyco_effect_unique_df, compare_unique_glyco_file(biomuta_df_chunk, glyco_effect_df)])
        acc_glyco_loss_unique_df = pd.concat([glyco_loss_unique_df, compare_unique_glyco_file(acc_chunk_df, glyco_loss_df)])
        acc_glyco_effect_unique_df = pd.concat([glyco_effect_unique_df, compare_unique_glyco_file(acc_chunk_df, glyco_effect_df)])

        

        
        
        # Add mutations for the specified accession
        accession_df = pd.concat([accession_df, acc_chunk_df])
    

    # Remove duplicate lines
    glyco_loss_overlap_df.drop_duplicates(keep='first', inplace=True)  
    glyco_effect_overlap_df.drop_duplicates(keep='first', inplace=True)
    glyco_loss_unique_df.drop_duplicates(keep='first', inplace=True)  
    glyco_effect_unique_df.drop_duplicates(keep='first', inplace=True)
    accession_df.drop_duplicates(keep='first', inplace=True) 
    

    # Finalize the glycosylation comparison
    glyco_report_dict = {}
    glyco_report_dict['Glycoloss shared'] = len(glyco_loss_overlap_df.index)
    glyco_report_dict['Glycoloss missing'] = len(glyco_loss_unique_df.index)
    glyco_report_dict['Glycoeffect shared'] = len(glyco_effect_overlap_df.index)
    glyco_report_dict['Glycoeffect missing'] = len(glyco_effect_unique_df.index)
    
    print(glyco_report_dict)
    
    # Export the shared and unique mutations for the biomuta vs. glyco datasets comparison
    glyco_loss_shared_path = output_folder + 'glyco_loss_overlap.csv'
    glyco_loss_unique_path = output_folder + 'glyco_loss_unique.csv'
    print("Exporting overall glycoloss shared and  mutations to: ")
    print(glyco_loss_shared_path)
    print(glyco_loss_unique_path)
    glyco_loss_overlap_df.to_csv(glyco_loss_shared_path, index = False)
    glyco_loss_unique_df.to_csv(glyco_loss_unique_path, index = False)

    glyco_effect_shared_path = output_folder + 'glyco_effect_overlap.csv'
    glyco_effect_unique_path = output_folder + 'glyco_effect_unique.csv'
    print("Exporting overall glycoeffect shared and missing mutations to:")
    print(glyco_effect_shared_path)
    print(glyco_effect_unique_path)
    glyco_effect_overlap_df.to_csv(glyco_effect_shared_path, index = False)
    glyco_effect_unique_df.to_csv(glyco_effect_unique_path, index = False)

    acc_glyco_loss_shared_path = output_folder + accession + '_glyco_loss_overlap.csv'
    acc_glyco_loss_unique_path = output_folder + accession + '_glyco_loss_unique.csv'
    print("Exporting glycoloss overlapping mutations for " + accession + " to: ")
    print(acc_glyco_loss_shared_path)
    print(acc_glyco_loss_unique_path)
    acc_glyco_loss_overlap_df.to_csv(acc_glyco_loss_shared_path, index = False)
    acc_glyco_loss_unique_df.to_csv(acc_glyco_loss_unique_path, index = False)

    acc_glyco_effect_shared_path = output_folder + accession + '_glyco_effect_overlap.csv'
    acc_glyco_effect_unique_path = output_folder + accession + '_glyco_effect_unique.csv'
    print("Exporting glycoeffect overlapping mutations for " + accession + " to:")
    print(acc_glyco_effect_shared_path)
    print(acc_glyco_effect_unique_path)
    acc_glyco_effect_overlap_df.to_csv(acc_glyco_effect_shared_path, index = False)
    acc_glyco_effect_unique_df.to_csv(acc_glyco_effect_unique_path, index = False)

    final_file_path = output_folder + accession + '.csv'
    print("Exporting all mutations for " + accession + " to " + final_file_path)
    accession_df.to_csv(final_file_path, index = False)



#def assess_glycosylation(ref_aa, alt_aa):
#    glycosylation_aa_list = ['T', 'S']

# Compare the AA mutations in biomuta to the glycosylation effect datasets from Glygen
def compare_shared_glyco_file(biomuta_df, glyco_df):

    # The compare field list is used to find overlapping data with the glyco datasets
    compare_field_list = [
        'uniprotkb_canonical_ac', 
        'aa_pos', 
        'ref_aa', 
        'alt_aa', 
        'do_name'
        ]
    
    # The output fields for the overlapping rows 
    output_field_list = [
        'uniprotkb_canonical_ac', 
        'aa_pos', 
        'ref_aa', 
        'alt_aa', 
        'do_name',
        'source'
        ]

    glyco_compare = glyco_df[compare_field_list].copy()
    
    biomuta_compare = biomuta_df[output_field_list].copy()
    
    shared_df = glyco_compare.merge(biomuta_compare,left_on=compare_field_list,right_on=compare_field_list,how='inner', indicator=True)

    return shared_df

def compare_unique_glyco_file(biomuta_df, glyco_df):

    # The compare field list is used to find overlapping data with the glyco datasets
    compare_field_list = [
        'uniprotkb_canonical_ac', 
        'aa_pos', 
        'ref_aa', 
        'alt_aa', 
        'do_name'
        ]
    
    # The output fields for the overlapping rows 
    output_field_list = [
        'uniprotkb_canonical_ac', 
        'aa_pos', 
        'ref_aa', 
        'alt_aa', 
        'do_name',
        ]

    glyco_compare = glyco_df[compare_field_list].copy()
    
    biomuta_compare = biomuta_df[output_field_list].copy()
    
    unique_df = glyco_compare.merge(biomuta_compare,left_on=compare_field_list,right_on=compare_field_list,how='left', indicator=True)
    unique_df.query("_merge == 'left_only'", inplace=True)

    return unique_df

def remove_duplicates(df):
    # Remove duplicate lines
    with_dup_total = len(df.index)
    df_out = df.drop_duplicates(keep='first')
    dup_number = len(df_out.index) - with_dup_total
    print("Removed " + str(dup_number) + " duplicate rows") 

                        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for csv combine script')
    parser.add_argument('--biomuta_csv', '-b',
                        help='An absolute path to the final biomuta mutations file')
    parser.add_argument('--accession', '-a',
                        help='The name of the uniprot accession to take mutations for', nargs='?', default='P00533-1')
    parser.add_argument('--glyco_file_folder', '-g',
                        help='An absolute path to the folder containing the glycosylation datasets from Glygen')
    parser.add_argument('--output_folder', '-o',
                        help='An absolute path to the final biomuta mutations file')
    args = parser.parse_args()

    main(args.biomuta_csv, args.accession, args.glyco_file_folder, args.output_folder)

#python accession_qc_check.py -b /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled/biomuta_v5.csv -g /mnt/c/Users/caule/github_general/biomuta/pipeline/qc_step4/glycosylation_datasets/ -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled/qc/