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
	

	uniprot_ann_file =  config_json["annotator"]["outdir"] + "/uniprot-ann.csv";
	seen = {"sitelist":{}}       
	with open(uniprot_ann_file, 'r') as FR:
                reader = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in reader:
			canon = row[0]
			if row[1] == "Glycosylation_Annotation":
				for pos in xrange(int(row[2]), int(row[3]) + 1): 
					if canon not in seen["sitelist"]:
						seen["sitelist"][canon] = []
					seen["sitelist"][canon].append(pos)
					
	call_dict = { True:{"+++":"no_change", "---":"loss"}, False:{"+++":"gain", "---":"no_change"}}
      
	out_file = config_json["annotator"]["outdir"] + "/netnglyc.2.csv"
	FW = open(out_file, "w") 
	chr_list = map(str,range(1,23)) + ["X","Y"]
        for chr_id in chr_list:
		in_file =  config_json["annotator"]["outdir"] + "/netnglyc.chr"+chr_id+".1.csv"
		with open(in_file, 'r') as FR:
                	reader = csv.reader(FR, delimiter=',', quotechar='"')
                	row_count = 0
                	for row in reader:
                        	row_count += 1
				canon = row[0]
				pos,refaa,altaa = row[1].split(":")
				r, motif = row[2].split(":")
				start_pos, end_pos = int(r.split("-")[0]), int(r.split("-")[1])
				label = row[-1]
				cond_list = [altaa != "*", label in ["+++", "---"], canon in seen["sitelist"]]

				if False not in cond_list:
					ann_flag = start_pos in seen["sitelist"][canon]
					call = call_dict[ann_flag][label]
					FW.write("%s,%s,%s,%s,%s,%s,%s\n" % (canon,pos,refaa,altaa,ann_flag,label,call))
				
	FW.close()
	

	



if __name__ == '__main__':
        main()








