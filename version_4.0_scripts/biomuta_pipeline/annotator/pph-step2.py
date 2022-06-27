import os,sys
import string
import csv
import json
from optparse import Option_parser
from Bio import Seq_iO
from Bio.Seq import Seq
import commands
import subprocess


__version__="1.0"
__status__ = "Dev"


#################################
def report_progress(progress_file, msg, open_new=False):

	msg = msg.strip()
        if open_new == True:
                with open(progress_file, "w") as FW:
                        FW.write(msg + "\n")
        else:
                with open(progress_file, "a") as FA:
                        FA.write(msg + "\n")
        return




###############################
def main():

	
	usage = "\n%prog  [options]"
	parser = Option_parser(usage,version="%prog " + __version__)
	parser.add_option("-i","--configfile",action="store",dest="configfile",help="Config file")



	(options,args) = parser.parse_args()
	for file in ([options.configfile]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	config_json = json.loads(open(options.configfile, "r").read())


	progress_file = config_json["annotator"]["outdir"] + "/pph-progress-step3.txt"


	seen = {}
        with open(config_json["annotator"]["outdir"] + "/precomputed.csv", "r") as FR:
                data_grid = csv.reader(FR, delimiter=',', quotechar='"')
                for row in data_grid:
			ac = row[0].split(".")[0]
			seen[ac] = True

        chr_list = map(str,range(1,23)) + ["X","Y"]
        for chr_id in chr_list:
                subs_file =  config_json["annotator"]["outdir"] + "/pph-subs.chr"+chr_id+".txt"
                out_file = config_json["annotator"]["outdir"] + "/pph-fsubs.chr"+chr_id+".txt"
		FW = open(out_file, "w")
		with open(subs_file, 'r') as FR:
                        reader = csv.reader(FR, delimiter='\t', quotechar='"')
                        row_count = 0
                        for row in reader:
				ac = row[0].strip()
				if ac in seen:
					FW.write("%s\n" % ("\t".join(row)))
		FW.close()
	
			


if __name__ == '__main__':
        main()








