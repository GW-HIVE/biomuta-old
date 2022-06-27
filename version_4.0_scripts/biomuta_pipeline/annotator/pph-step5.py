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


		
	uniprot_file = config_json["uniprotfile"]
	pph_uniprot_file = config_json["pphuniprotfile"]



 	seen = {"canon":{}}

        with open(config_json["idmapping"], "r") as FR:
                data_grid = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in data_grid:
                        isoform = row[0]
                        seen["canon"][isoform] = True



	mapped2genome = {}
        with open(config_json["pepid2isoformac"], "r") as FR:
                data_grid = csv.reader(FR, delimiter=',', quotechar='"')
                for row in data_grid:
                        pep_id, isoform_list = row[0], row[2].split("|")
                        if row[1] != "ignore":
                                for isoform in isoform_list:
                                        mapped2genome[isoform] = pep_id

	seq2aclist = {}
        for record in Seq_iO.parse(config_json["pphuniprotfile"], "fasta"):
		ac = record.id.split("|")[1]
                seq = str(record.seq.upper())
                if seq not in seq2aclist:
                        seq2aclist[seq] = {"list1":[], "list2":[]}
                seq2aclist[seq]["list1"].append(ac)

        for record in Seq_iO.parse(config_json["uniprotfile"], "fasta"):
                isoform = record.id.split("|")[1]
		seq = str(record.seq.upper())
		if seq not in seq2aclist:
			seq2aclist[seq] = {"list1":[], "list2":[]}
		if isoform in mapped2genome:
			seq2aclist[seq]["list2"].append(isoform)

	ac2isoformlist = {}
	for seq in seq2aclist:
		list1, list2 = seq2aclist[seq]["list1"], seq2aclist[seq]["list2"]
		for ac in list1:
			ac2isoformlist[ac] = list2

	

	out_file = config_json["annotator"]["outdir"] + "/pph-final.csv"
	FW = open(out_file, "w")
	chr_list = map(str,range(1,23)) + ["X","Y"]
	for chr_id in chr_list:
		predictions_file =  config_json["annotator"]["outdir"] + "/pph-predictions.chr"+chr_id+".txt"
		with open(predictions_file, 'r') as FR:
               		reader = csv.reader(FR, delimiter='\t', quotechar='"')
               		row_count = 0
               		for row in reader:
           			row_count += 1
				if row_count == 1:
					continue
				ac,pos,ref,alt = row[0].strip(),row[1].strip(),row[2].strip(),row[3].strip()
				pred,prob = row[11].strip(),row[15].strip()
				if ac in ac2isoformlist:
					canon = ""
					for isoform in ac2isoformlist[ac]:
						if isoform in seen["canon"]:
							canon = isoform
							break
					if canon != "":
						FW.write("%s,%s,%s,%s,%s,%s\n" % (canon,pos,ref,alt, pred,prob))
	FW.close()







if __name__ == '__main__':
        main()








