import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import csv


__status__ = "Dev"
__version__ = "4.0"






##############################################
def main():

        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-i", "--configfile", action = "store", dest = "configfile", help = "Config file")
	parser.add_option("-s", "--src_type", action = "store", dest = "src_type", help = "src_type") 


	(options,args) = parser.parse_args()
        for file in ([options.configfile, options.src_type]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

	json_file = options.configfile
        config_json = json.loads(open(json_file).read())
	src_type = options.src_type

	#chr_list = map(str,range(1,23)) + ["X","Y"]
        chr_list = ["14"]

	input_dir = config_json["datareporter"][src_type]["inputdir"]
	if src_type == "tcga":
		pattern = input_dir + "/*/*/*.vcf"
	else:
		pattern = input_dir + "/*.vcf"
	file_list = glob.glob(pattern)
	for chr_id_out in chr_list:
		file_count = 0
		mutid_dict = {"cosmic":{}, "tcga":{}, "icgc":{}}
		print "started %s of %s ... " % (chr_id_out, src_type) 
		for vcf_file in file_list:
			patient_id = vcf_file.split("/")[-1].split(".")[0]
			file_count += 1
			with open(vcf_file, 'r') as FR:
               			row_count = 0
				for line in FR:
					row_count += 1
					row = line.strip().split("\t")
					if row[0][0:1] == "#":	
						continue
					chr_id = row[0]
					if chr_id not in chr_list:
						continue
					if src_type == "tcga" and chr_id != chr_id_out:
						continue                                        
       					combo_id = "%s,%s,%s,%s,%s" % (row[0],row[1],row[3],row[4],patient_id)
                                        if row[0] == "14" and row[1] == "58606016":
                                            print combo_id, vcf_file




if __name__ == '__main__':
	main()


