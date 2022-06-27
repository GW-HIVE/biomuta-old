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

	

	chr_list = map(str,range(1,23)) + ["X","Y"]

	for chr_id in chr_list:
		features_file =  config_json["annotator"]["outdir"] + "/pph-features-step2.chr"+chr_id+".txt"
		predictions_file =  config_json["annotator"]["outdir"] + "/pph-predictions.chr"+chr_id+".txt"
		log_file = config_json["annotator"]["outdir"] + "/pph-step3.chr"+chr_id+".log"
		cmd = config_json["binaries"]["runweka"] + " " + features_file 
		cmd += " 1>" + predictions_file + " 2>" + log_file   
		#print cmd
		commands.getoutput(cmd)




if __name__ == '__main__':
        main()








