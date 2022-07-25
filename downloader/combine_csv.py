'''
Input:
########
    * -i : The folder containing csv mutation files to combine

Output:
########
    * A csv file combining all csv files in a given folder

Usage:
########
    * python combine_csv.py -h

    *Gives a description of the neccessary commands

    * python combine_csv.py -i <path/>

    *Runs the script with the given folder and combines all csv files in that folder

'''

import argparse
import os
import glob
import pandas as pd



def main(input_folder):
    
    # Grab all of the csv files in the specified folder
    final_mutation_files = os.path.join(input_folder, "*.csv")
    final_mutation_files = glob.glob(final_mutation_files)

    final_df = pd.DataFrame()
    
    for file in final_mutation_files:
        print("Adding " + file)
        temp_df = pd.read_csv(file, dtype=str)
        final_df = pd.concat([final_df, temp_df])


    # Combine the files into a master dataframe
    #final_df = pd.concat(map(pd.read_csv, final_mutation_files), ignore_index=True)
    
    # Convert the 
    #final_df['aa_pos']

    final_file_path = input_folder + "/biomuta_v5.csv"
    print("Exporting mapped file to " + final_file_path)
    final_df.to_csv(final_file_path, index = False)

                        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for csv combine script')
    parser.add_argument('--input_folder', '-i',
                        help='An absolute path to folder containing csv files to combine')
    args = parser.parse_args()

    main(args.input_folder)

#python combine_csv.py -i /mnt/c/Users/caule/OncoMX/biomuta/v-5.0/compiled/