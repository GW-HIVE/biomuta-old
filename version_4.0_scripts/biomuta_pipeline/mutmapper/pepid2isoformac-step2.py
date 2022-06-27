import os,sys
import string
import csv
import json
from optparse import Option_parser
from Bio import Seq_iO
from Bio.Seq import Seq
from Bio import pairwise2


__version__="1.0"
__status__ = "Dev"


#################################
def report_progress(progress_file, msg, open_new=False):

        if open_new == True:
                with open(progress_file, "w") as FW:
                        FW.write(msg)
        else:
                with open(progress_file, "a") as FA:
                        FA.write(msg)
        return



###############################
def main():

	usage = "\n%prog  [options]"
	parser = Option_parser(usage,version="%prog " + __version__)
	parser.add_option("-i","--configfile",action="store",dest="configfile",help="NT file")


	(options,args) = parser.parse_args()
	for file in ([options.configfile]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	config_json = json.loads(open(options.configfile, "r").read())
	pep_file = config_json["pepfile"]
	uniprot_file = config_json["uniprotfile"]
	out_file = config_json["mutmapper"]["outdir"] + "/map2uniprot.ranges.csv"        

 
	#Load cds sequences
	ens_pep_seq_hash = {}
	for record in Seq_iO.parse(pep_file, "fasta"):
		ens_pep_seq_hash[record.id] = record.seq.upper()
        

        uniprot_pep_seq_hash = {}
        for record in Seq_iO.parse(uniprot_file, "fasta"):
                isoform_ac = record.id.split("|")[1]
		uniprot_pep_seq_hash[isoform_ac] = record.seq.upper()


	FW = open(out_file, "w")
	FW.write("peptide_id,isoform_ac,ranges\n")
        with open(config_json["pepid2isoformac"], "r") as FR:
                data_grid = csv.reader(FR, delimiter=',', quotechar='"')
                row_count = 0
                for row in data_grid:
			pep_id, ac_list = row[0], row[2].split("|")
			if row[1] in ["exact", "ignore"]:
				continue
			for isoform_ac in ac_list:
				if isoform_ac not in uniprot_pep_seq_hash:
					continue
				seq1, seq2 = ens_pep_seq_hash[pep_id], uniprot_pep_seq_hash[isoform_ac]
				aln_list = pairwise2.align.globalxx(seq1, seq2)
				pos1, pos2 = 0, 0
				prev_pos1, prev_pos2 = 0, 0
				prev_aa1, prev_aa2 = "-", "-"
				start_pos1, start_pos2, end_pos1, end_pos2 = 0, 0, 0, 0
				range_list = []
				alnlen = 0
				for i in xrange(0, len(aln_list[0][0])):
					aa1, aa2 = aln_list[0][0][i], aln_list[0][1][i]
					pos1 += 1 if aa1 != "-" else 0
					pos2 += 1 if aa2 != "-" else 0
					c1 = aa1 != "-" and aa2 != "-"
                                	c2 = prev_aa1 == "-" or prev_aa2 == "-"
                                	if c1 and c2:
                                        	start_pos1, start_pos2 = pos1, pos2
					c1 = aa1 == "-" or aa2 == "-"
					c2 = prev_aa1 != "-" and prev_aa2 != "-"
					if (c1 and c2) or pos1 == len(seq1):
						end_pos1, end_pos2 = prev_pos1, prev_pos2
						alnlen += end_pos1 - start_pos1 + 1
						r = "%s-%s:%s-%s"%(start_pos1,end_pos1,start_pos2,end_pos2)
						range_list.append(r)
					prev_aa1 = aa1
					prev_aa2 = aa2
					prev_pos1 = pos1
					prev_pos2 = pos2
				ratio1 = float(alnlen)/len(seq1) 
				ratio2 = float(alnlen)/len(seq2)
				if ratio1 > 0.7 and len(range_list) <= 10:
					FW.write("%s,%s,%s\n" % (pep_id, isoform_ac, "|".join(range_list)))
	FW.close()







if __name__ == '__main__':
        main()








