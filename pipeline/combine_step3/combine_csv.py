'''
Input:
########
    * -i : The folder containing csv mutation files to combine
    * -o : The folder to output the combined mutation file

Output:
########
    * A csv file combining all csv files in a given folder

Usage:
########
    * python combine_csv.py -h

    *Gives a description of the neccessary commands

    * python combine_csv.py -i <path/> -o <path/>

    *Runs the script with the given folder and combines all csv files in that folder

'''

import argparse
import os
import glob
import pandas as pd
import re



def get_next_version_number(output_folder, base_filename):
    """
    Finds the highest version number in the output folder for files matching the base filename and returns the next version number.
    """
    # Regex pattern to match "biomuta_vX.csv"
    pattern = re.compile(rf"{re.escape(base_filename)}_v(\d+)\.csv$")

    max_version = 0

    # Iterate through files in the output folder
    for filename in os.listdir(output_folder):
        match = pattern.match(filename)
        if match:
            # Extract the version number from the filename
            version = int(match.group(1))
            if version > max_version:
                max_version = version

    # Return the next version number
    return max_version + 1


def main(input_folder, output_folder):
    
    # Grab and combine all of the csv files in the specified folder
    final_mutation_files = os.path.join(input_folder, "*.csv")
    final_mutation_files = glob.glob(final_mutation_files)

    final_df = pd.DataFrame()
    
    for file in final_mutation_files:
        print("Adding " + file)
        temp_df = pd.read_csv(file, dtype=str)
        final_df = pd.concat([final_df, temp_df])

   
    # Remove duplicate lines
    with_dup_total = len(final_df)
    final_df.drop_duplicates(keep='first', inplace=True)
    dup_number = len(final_df.index) - with_dup_total
    if dup_number > 0:
        print("Removed " + str(dup_number) + " duplicate rows") 

    base_filename = "biomuta"
    next_version = get_next_version_number(output_folder, base_filename)
    new_filename = f"{base_filename}_v{next_version}.csv"
    final_file_path = os.path.join(output_folder, new_filename)
    print("Exporting mapped file to " + final_file_path)
    final_df.to_csv(final_file_path, index = False, quoting=1)

                        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Commands for csv combine script')
    parser.add_argument('--input_folder', '-i',
                        help='An absolute path to a folder containing csv files to combine')
    parser.add_argument('--output_folder', '-o',
                        help='An absolute path to a folder where the final combined mutation file will be written')
    args = parser.parse_args()

    main(args.input_folder, args.output_folder)

#python3 combine_csv.py -i /data/shared/repos/biomuta-old/generated_datasets/compiled/source_mutation_files -o /data/shared/repos/biomuta-old/generated_datasets/compiled