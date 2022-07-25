'''
Input:
########
    * -b : The final biomuta mutation file to extract mutations from
    * -a : The file containing a list of uniprot accessions for mutation extraction
    * -o : The output folder to export mutations for selected accessions

Output:
########
    * A csv file per gene, containing mutations extracted from the main biomuta mutation file

Usage:
########
    * python combine_csv.py -h

    *Gives a description of the neccessary commands

    * python combine_csv.py -b <path/final_biomuta_file.csv> -a <path/gene_list.csv> -o <path/>

    *Runs the script with the final biomuta mutation file and extracts mutations per accession given in the accession list

'''

import argparse
import pandas as pd
import csv



def main(biomuta_csv, accession_list, output_folder):

    with open(accession_list, "r") as accession_list_handle:
        accessions = csv.reader(accession_list_handle)

        # Save the list of accessions
        acc_list = []
        for acc in accessions:
            acc_list.append(acc)
    
    # Set up the chunk size for better memory usage and script speed
    biomuta_df_iterator = pd.read_csv(biomuta_csv, dtype=str, chunksize=1000000)
    
    # Set up the gene specific dataframes
    for acc in acc_list:

        acc_df_name = acc + "_df"
        pd.DataFrame(acc_df_name)
        

    for i, biomuta_df_chunk in enumerate(biomuta_df_iterator):

        print("Processing chunk number " + str(i+1))

        # Add rows containing the gene to the gene specific dataframe
        for acc in acc_list:

            acc_df_name = acc + "_df"

            temp_df = biomuta_df_chunk[biomuta_df_chunk['uniprot_canonical_ac'] == acc]
    
            acc_df_name = pd.concat[acc_df_name, temp_df]

    for acc in acc_list:

        acc_df_name = acc + "_df"
    
        final_file_path = output_folder + acc + '.csv'
        print("Exporting mutations for " + acc + " to " + final_file_path)
        acc_df_name.to_csv(final_file_path, index = False)

                        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for csv combine script')
    parser.add_argument('--biomuta_csv', '-b',
                        help='An absolute path to the final biomuta mutations file')
    parser.add_argument('--accession_list', '-a',
                        help='An absolute path to the file containing a list of uniprot accessions')
    parser.add_argument('--biomuta_csv', '-b',
                        help='An absolute path to the final biomuta mutations file')
    args = parser.parse_args()

    main(args.biomuta_csv, args.qc_genes_list, args.output_folder)

#python combine_csv.py -i /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled/