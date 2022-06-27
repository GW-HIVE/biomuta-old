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


	len_hash = {}
	seq2aclist = {}
	pepid2trsid = {}	
        for record in Seq_iO.parse(config_json["pepfile"], "fasta"):
                trs_id = record.description.split(" ")[4].split(":")[1]
		pep_id = record.id
		pepid2trsid[pep_id] = trs_id
                seq = str(record.seq.upper())
                if seq not in seq2aclist:
                        seq2aclist[seq] = {"list1":[], "list2":[]}
                seq2aclist[seq]["list1"].append(pep_id)
		len_hash[pep_id] = len(seq)

	ac2isoformlist = {}
        for record in Seq_iO.parse(config_json["uniprotfile"], "fasta"):
                isoform = record.id.split("|")[1]
		seq = str(record.seq.upper())
		if seq not in seq2aclist:
			seq2aclist[seq] = {"list1":[], "list2":[]}
		seq2aclist[seq]["list2"].append(isoform)
		ac = isoform.split("-")[0]
		if ac not in ac2isoformlist:
			ac2isoformlist[ac] = []
		ac2isoformlist[ac].append(isoform)		



	pepid2isoformac = {}
	for seq in seq2aclist:
		c1 = len(seq2aclist[seq]["list1"]) > 0
		c2 = len(seq2aclist[seq]["list2"]) > 0
		if c1 and c2:
			for pep_id in seq2aclist[seq]["list1"]:
                                pepid2isoformac[pep_id] = "exact,"+"|".join(seq2aclist[seq]["list2"])



        cmd = "mysql -u anonymous -h ensembldb.ensembl.org homo_sapiens_core_75_37 "
        cmd += "--execute \"SELECT A.stable_id, C.dbprimary_acc FROM translation A, "
        cmd += "object_xref B, xref C, external_db D WHERE A.translation_id = B.ensembl_id "
        cmd += "AND B.ensembl_object_type = 'Translation' AND B.xref_id = C.xref_id AND "
        cmd += "C.external_db_id = D.external_db_id AND D.db_name = 'Uniprot/SPTREMBL'\""
        lines = commands.getoutput(cmd).split("\n")

        extra_list = {}
        for i in xrange(1, len(lines)):
                pep_id, ac = lines[i].strip().split("\t")
		if pep_id not in pepid2isoformac:
                       	if ac in ac2isoformlist:
				if pep_id not in extra_list:
					extra_list[pep_id] = []
				for isoform in ac2isoformlist[ac]:
					extra_list[pep_id].append(isoform)

        for pep_id in extra_list:
                flag = "nonexact" if len_hash[pep_id] < 5000 else "ignore"
		pepid2isoformac[pep_id] = flag + "," + "|".join(extra_list[pep_id])


	FW = open(config_json["pepid2isoformac"], "w")
        for pep_id in pepid2isoformac:
                FW.write("%s,%s\n" % (pep_id,  pepid2isoformac[pep_id]))
	FW.close()






	
	
			


if __name__ == '__main__':
        main()








