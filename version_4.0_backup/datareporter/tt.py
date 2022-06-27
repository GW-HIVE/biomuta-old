import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import csv


__status__ = "Dev"
__version__ = "4.0"






#######################################
def load_do_slim (doid_file, doslimid_hash, doslimname_hash, exclude_list):


        with open(doid_file, 'r') as FR:
                csv_grid = csv.reader(FR, delimiter='\t', quotechar='"')
                for row in csv_grid:
                        doslim_id = row[0].split(" ")[0].split(":")[1]
                        if doslim_id in exclude_list:
                                continue
                        doslimname_hash[doslim_id] = row[0]
                        doslimid_hash[row[0]] = doslim_id
                        doslimid_hash[doslim_id] = doslim_id
			for val in row[1].split("|"):
                                doslimid_hash[val] = doslim_id
                        for j in [2,3]:
                                val = row[j].strip()
                                if val != "":
                                        doslimid_hash[val] = doslim_id
        return



##############################################
def main():

        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-i", "--configfile", action = "store", dest = "configfile", help = "Config json file path")
        
	(options,args) = parser.parse_args()
        for file in ([options.configfile]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

	src_type = "tcga"
	json_file = options.configfile
        config_json = json.loads(open(json_file).read())

	input_dir = config_json["datareporter"][src_type]["inputdir"]
	output_dir = config_json["datareporter"][src_type]["outputdir"]

	global count_hash
	global seen_hash
	global cosmicid2doid

		
	count_hash = {
		"badrows":{},
		"variants_badrows":0,
		"variants_total":0,
		"variants_valid":0,
		"variants_valid_unique":0
	}
	doslimid_hash = {}
	doslimname_hash = {}



	#Load bad rows
	bad_row_hash = {}	
	bad_rows_file = output_dir + "/badrows.txt"
	with open(bad_rows_file, 'r') as FR:
		csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
		for row in csv_grid:
			if row[0] not in count_hash["badrows"]:
				count_hash["badrows"][row[0]] = 0	
			count_hash["badrows"][row[0]] += 1
			if row[2] not in bad_row_hash:
				bad_row_hash[row[2]] = {}
			bad_row_hash[row[2]][row[1]] = True

	#Load DOID mappings
	doid_file = config_json["datareporter"]["doid"]["outputdir"] + "/doid-mapping.reduced.tsv"
	load_do_slim (doid_file, doslimid_hash, doslimname_hash, [])
        #Repeat loading to overwrite mapping to 305 if more specific mapping exists
	load_do_slim (doid_file, doslimid_hash, doslimname_hash, ["305"])




	chr_list1 = map(str,range(1,23)) + ["X","Y"]
	chr_list2 = ["1"]
	if src_type == "tcga":
		pattern = input_dir + "/*/*/*.vcf"
	else:
		pattern = input_dir + "/*.vcf"
	file_list = glob.glob(pattern)



	nx = 0
	cancer_type_hash = {}
	chr_list = chr_list2 if src_type == "tcga" else chr_list1
	for chr_id in chr_list:
		freq_hash = {}
        	seen_hash = {"comboid":{}, "nondoslim":{}}
		newdoslim_id = 0
		file_count = 0
		for vcf_file in file_list:
			cancer_type = vcf_file.split("/")[-3]
			file_name = vcf_file.split("/")[-1]
			file_count += 1
			with open(vcf_file, 'r') as FR:
                		row_count = 0
				for line in FR:
					row_count += 1
					row = line.strip().split("\t")
					if row[0][0:1] == "#":	
						continue
					if src_type!= "tcga" and row[0] != chr_id:
						continue
					chr_id = row[0] if src_type == "tcga" else chr_id 
					if chr_id not in chr_list1:
						continue                                        
					count_hash["variants_total"] += 1
					if vcf_file in bad_row_hash:
						if str(row_count) in bad_row_hash[vcf_file]:
							count_hash["variants_badrows"] += 1
							continue 
					pf = True
        				do_id = ""
        				print row



if __name__ == '__main__':
	main()


