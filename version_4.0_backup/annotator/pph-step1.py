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
	parser.add_option("-c","--chrid",action="store",dest="chrid",help="Chr ID (1-22,X,Y,MT) ")



	(options,args) = parser.parse_args()
	for file in ([options.configfile, options.chrid]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	chr_id = options.chrid
	config_json = json.loads(open(options.configfile, "r").read())
	merged_file = config_json["mutmapper"]["outdir"] + "/mutlist-merged.csv"


	progress_file = config_json["annotator"]["outdir"] + "/pph-progress-step1."+chr_id+".txt"
	snps_file =  config_json["annotator"]["outdir"] + "/pph-snps.chr"+chr_id+".txt"
	FW = open(snps_file, "w")
	start_flag = False
       	with open(merged_file, 'r') as FR:
               	reader = csv.reader(FR, delimiter=',', quotechar='"')
               	row_count = 0
               	for row in reader:
                       	row_count += 1
                       	if row_count == 1:
                               	continue
			if row[1] != chr_id:
				if start_flag == True:
					break
				else:
					continue
			if row[-1] == "noncoding":
				continue
			start_fag = True
			chr_lbl = "M" if row[1] == "MT" else row[1]
			FW.write("chr%s:%s %s/%s\n" % (chr_lbl,row[2],row[3],row[4]))
			if row_count%1000000 == 0:
				report_progress(progress_file, "Parsed %s rows of chr%s" % (row_count, chr_id))
	FW.close()


	subs_file =  config_json["annotator"]["outdir"] + "/pph-subs.chr"+chr_id+".txt"
	features_file =  config_json["annotator"]["outdir"] + "/pph-features-step1.chr"+chr_id+".txt"
	log_file = config_json["annotator"]["outdir"] + "/pph-step1.chr"+chr_id+".log"
	cmd = config_json["binaries"]["mapsnps"] + " -g hg19 -m -U -y " + subs_file + " " + snps_file
	cmd += " 1>" + features_file + " 2>" + log_file   
	#print cmd
	commands.getoutput(cmd)




if __name__ == '__main__':
        main()








