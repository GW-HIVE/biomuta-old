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

    print(glyco_file_list)
    
    # Set up dataframes to hold the glycosylation comparison
    for file in glyco_file_list:
        if re.search(r'loss', file):
            print('Glycoloss file: ' + file)
            glyco_loss_df = pd.read_csv(file, dtype = str)
        elif re.search(r'effect', file):
            print('Glycoeffect file ' + file)
            glyco_effect_df = pd.read_csv(file, dtype=str)
        else:
            print('Glycosylation dataset ' + file + ' not used for comparison')
    
    glyco_loss_overlap_df = pd.DataFrame()
    glyco_effect_overlap_df = pd.DataFrame()
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
        glyco_loss_overlap_df = pd.concat([glyco_loss_overlap_df, compare_glyco_file(biomuta_df_chunk, glyco_loss_df)])
        glyco_effect_overlap_df = pd.concat([glyco_effect_overlap_df, compare_glyco_file(biomuta_df_chunk, glyco_effect_df)])
        acc_glyco_loss_overlap_df = pd.concat([glyco_loss_overlap_df, compare_glyco_file(acc_chunk_df, glyco_loss_df)])
        acc_glyco_effect_overlap_df = pd.concat([glyco_effect_overlap_df, compare_glyco_file(acc_chunk_df, glyco_effect_df)])

        # Add mutations for the specified accession
        accession_df = pd.concat([accession_df, acc_chunk_df])
    

    # Finalize the glycosylation comparison
    glyco_report_dict = {}
    glyco_report_dict['Glycoloss shared'] = len(glyco_loss_overlap_df.index)
    glyco_report_dict['Glycoloss missing'] = len(glyco_loss_df.index) - len(glyco_loss_overlap_df.index)
    glyco_report_dict['Glycoeffect shared'] = len(glyco_effect_overlap_df.index)
    glyco_report_dict['Glycoeffect missing'] = len(glyco_effect_df.index) - len(glyco_effect_overlap_df.index)

    glyco_loss_path = output_folder + accession + '_glyco_loss_overlap.csv'
    print("Exporting glycoloss overlapping mutations for to " + glyco_loss_path)
    glyco_loss_overlap_df.to_csv(glyco_loss_path, index = False)

    glyco_effect_path = output_folder + accession + '_glyco_effect_overlap.csv'
    print("Exporting glycoeffect overlapping mutations for to " + glyco_loss_path)
    glyco_effect_overlap_df.to_csv(glyco_effect_path, index = False)
        
    final_file_path = output_folder + accession + '.csv'
    print("Exporting mutations for " + accession + " to " + final_file_path)
    accession_df.to_csv(final_file_path, index = False)



#def assess_glycosylation(ref_aa, alt_aa):
#    glycosylation_aa_list = ['T', 'S']

# Compare the AA mutations in biomuta to the glycosylation effect datasets from Glygen
def compare_glyco_file(biomuta_df, glyco_df):
    compare_field_list = [
        'uniprotkb_canonical_ac', 
        'aa_pos', 
        'ref_aa', 
        'alt_aa', 
        'do_name'
        ]
    
    output_field_list = [
        'uniprotkb_canonical_ac', 
        'aa_pos', 
        'ref_aa', 
        'alt_aa', 
        'do_name'
        ]

    glyco_compare = glyco_df[compare_field_list].copy()
    
    biomuta_compare = biomuta_df[output_field_list].copy()
    
    shared_df = glyco_compare.merge(biomuta_compare,left_on=compare_field_list,right_on=compare_field_list,how='left')

    return shared_df
                        
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