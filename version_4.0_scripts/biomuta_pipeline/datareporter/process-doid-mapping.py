import os,sys
import string
from optparse import Option_parser
import glob
import json
import csv


__version__="4.0"
__status__ = "Dev"


########################################

def main():
	
	usage = "\n%prog  [options]"
        parser = Option_parser(usage,version="%prog " + __version__)
	parser.add_option("-i", "--configfile", action = "store", dest = "configfile", help = "Config json file path")
        
	(options,args) = parser.parse_args()
        for file in ([options.configfile]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

	
	config_file = options.configfile
	config_json = json.loads(open(config_file).read())
	in_file = config_json["datareporter"]["doid"]["input_dir"] + "/doid-mapping.tsv"
	out_file = config_json["datareporter"]["doid"]["output_dir"] + "/doid-mapping.reduced.tsv"

	FW = open(out_file, "w")
	
	with open(in_file, 'r') as FR:
               	csv_grid = csv.reader(FR, delimiter='\t', quotechar='"')
               	line_count = 0
		fieldname2index = {}
               	for row in csv_grid:
			line_count += 1
			if line_count == 1:
				for j in xrange(0, len(row)):
					fieldname2index[row[j]] = j
				continue
			new_row = []
			do_id_list = []
			for j in [3,7,9]:
				if row[j].strip() != "-":
					for x in row[j].split("DOID:")[1:]:
						do_id_list.append(x.split("/")[0].strip())
			short_name = ""
			if row[12] in ["TCGA", "ICGC", "TARGET", "ICGC-TARGET" , "ICGC-TCGA"]:
				short_name = row[11].split(" ")[-1].strip()[1:-1]
			
			val1 = row[1]
			val2 = "|".join(do_id_list)
			val3 = short_name
			val4 = row[11]
			if val1 != "" and (val2 != "" or val3 != ""):
				FW.write("%s\t%s\t%s\t%s\n" %  (val1, val2, val3, val4))


	FW.close()

if __name__ == '__main__':
        main()





