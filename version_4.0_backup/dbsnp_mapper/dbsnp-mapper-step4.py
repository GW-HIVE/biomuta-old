import os,sys
import string
import csv
import json
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq


__version__="1.0"
__status__ = "Dev"




###############################
def main():

        usage = "\n%prog  [options]"
        parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-i","--configfile",action="store",dest="configfile",help="NT file")


        (options,args) = parser.parse_args()
        for file in ([options.configfile]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

        config_json = json.loads(open(options.configfile, "r").read())

        dbsnp_in_file = "outdir/dbsnp.3.csv"
        dbsnp_out_file = "outdir/dbsnp.4.csv"

	FW = open(dbsnp_out_file, "w")
	with open(dbsnp_in_file, "r") as FR:
                data_grid = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in data_grid:
			row_count += 1
			if row[10] == "NSM" and row[-2] == row[-1]:
				continue
			if row[12] == "SYN" and row[-2] != row[-1]:
				continue
			FW.write("%s\n" % ("\t".join(row)))
    	FW.close()





if __name__ == '__main__':
        main()








