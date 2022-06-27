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
 

	(options,args) = parser.parse_args()
        for file in ([options.configfile]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

	json_file = options.configfile
        config_json = json.loads(open(json_file).read())


	chr_list = map(str,range(1,23)) + ["X","Y"]
        #chr_list = ["1"]

	pattern = config_json["mutmapper"]["outdir"] + "/mutid-tcga-chr*.csv"
	file_list = glob.glob(pattern)
	data_frame = {}
	src_list = []
	for in_file in file_list:
		file_name = in_file.split("/")[-1]
		src_type = file_name.split("-")[1]
		src_list.append(src_type)
		with open(in_file, 'r') as FR:
               		row_count = 0
			for line in FR:
				row_count += 1
				row = line.strip().split(",")
				id_list = row[-1].split("|")
				if len(id_list) > 1:
					for j in xrange(0, len(id_list)):
						id_list[j] = id_list[j].split("t")[0]
					id_list = sorted(set(id_list))
					if len(id_list) > 1:
						print  "%s,%s" % (",".join(row[:-1]), "|".join(id_list))





if __name__ == '__main__':
	main()


