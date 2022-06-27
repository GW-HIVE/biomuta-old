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

	out_file = config_json["annotator"]["outdir"] + "/netnglyc.chr"+chr_id+".1.csv"
        progress_file = config_json["annotator"]["outdir"] + "/netnglyc-step1-progress-chr"+chr_id+".txt"
        report_progress(progress_file, "Started netnglyc-step1 ...", True)
	uniprot_ann_file =  config_json["annotator"]["outdir"] + "/uniprot-ann.csv";


	seen = {"signalp":{}, "netnglyc":{} }       
	with open(uniprot_ann_file, 'r') as FR:
                reader = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in reader:
			if row[1] == "Signal_Peptide_Annotation":
				canon = row[0]
				seen["signalp"][canon] = True
	
	canon_seq_hash = {}
	for record in Seq_iO.parse(config_json["canonicaluniprotfile"], "fasta"):
                seq_id = record.id.split("|")[1]
		seq = str(record.seq.upper())
		canon_seq_hash[seq_id] = seq
	

	with open(out_file, "w") as FW:
		FW.write("")

	mutmap_file =  config_json["mutmapper"]["outdir"] + "/mutmap.2.csv"
        with open(mutmap_file, 'r') as FR:
                reader = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in reader:
                        row_count += 1
                        if row_count == 1:
                                continue
                        if row[1] != chr_id:
				continue
			mut_index,chrom,pos,ref,alt = row[0:5]
			refaa,altaa = row[12], row[13]
                        isoform_list,pos_in_uniprot_ac = row[14].split("|"), row[15]
			if refaa == altaa:
				continue
			if row[14] != "-":
				for isoform in isoform_list:
					if isoform in seen["signalp"]:
						key = "%s,%s,%s,%s" % (isoform,pos_in_uniprot_ac,refaa,altaa)
						if key in seen["netnglyc"]:
							continue
						seen["netnglyc"][key] = True
						pos = int(pos_in_uniprot_ac)
						variant_seq = canon_seq_hash[isoform][:pos-1] + altaa 
						variant_seq += canon_seq_hash[isoform][pos:]
						seq_file = config_json["annotator"]["outdir"] + "/netnglyc-batch-chr"+chr_id+".fasta"
						FW1 = open(seq_file, "w")
    						FW1.write(">%s\n%s\n" % (isoform,variant_seq))
						FW1.close()
						cmd = config_json["binaries"]["netnglyc"] + " " + seq_file
        					outlines = commands.getoutput(cmd).split("\n")
						if outlines[-2].find("---------------") == -1:
							continue
						i = -3
        					while outlines[i].find("---------------") == -1:
                					outrow = outlines[i].split()
							motif = outrow[2].replace("-", "")
							start_pos = int(outrow[1])
                					end_pos = start_pos + len(motif) - 1
							if pos >= start_pos and pos <= end_pos:
								r = "%s:%s:%s,%s-%s:%s" % (pos,refaa,altaa,start_pos, end_pos,motif)
								out_list = [outrow[0],r] + outrow[3:]
								with open(out_file, "a") as FA:
									FA.write("%s\n" % (",".join(out_list)))
							i -= 1
			if row_count%100000 == 0:
				report_progress(progress_file, "Processed mutation %s" % (row_count))

	



if __name__ == '__main__':
        main()








