'''
Input:
########
    * -b : The final biomuta mutation file to extract mutations from
    * -a : The name of the uniprot accession for mutation extraction
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



def main(biomuta_csv, accession, output_folder):

    # Set up a dataframe to hold on mutations for the accession
    accession_df = pd.DataFrame()
    
    # Set up the chunk size for better memory usage and script speed
    biomuta_df_iterator = pd.read_csv(biomuta_csv, dtype=str, chunksize=1000000)
    
    # Set up dictionary containing a list of mutations for each accession

    for i, biomuta_df_chunk in enumerate(biomuta_df_iterator):

        print("Processing chunk number " + str(i+1))

        chunk_df = biomuta_df_chunk[biomuta_df_chunk['uniprotkb_canonical_ac'] == accession]
    
        accession_df = pd.concat([accession_df, chunk_df])
    
    final_file_path = output_folder + accession + '.csv'
    print("Exporting mutations for " + accession + " to " + final_file_path)
    accession_df.to_csv(final_file_path, index = False)

                        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for csv combine script')
    parser.add_argument('--biomuta_csv', '-b',
                        help='An absolute path to the final biomuta mutations file')
    parser.add_argument('--accession', '-a',
                        help='The name of the uniprot accession to take mutations for')
    parser.add_argument('--output_folder', '-o',
                        help='An absolute path to the final biomuta mutations file')
    args = parser.parse_args()

    main(args.biomuta_csv, args.accession, args.output_folder)

#python accession_qc_check.py -b /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled/biomuta_v5.csv -a P00533-1 -o /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled/qc/