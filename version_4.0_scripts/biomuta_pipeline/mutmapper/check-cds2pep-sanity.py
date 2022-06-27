import os,sys
import string
import json
from optparse import Option_parser
from Bio import Seq_iO
from Bio.Seq import Seq


__version__="1.0"
__status__ = "Dev"


###############################
def my_translator(seq, codon_dic):

	protein_seq = ''
	for i in range(0,len(seq),3):
		if seq[i:i+3] in codon_dic:
			if i == len(seq) - 3 and codon_dic[seq[i:i+3]] == "*":
				protein_seq += ""
			elif codon_dic[seq[i:i+3]] == "*":
				protein_seq += "*"
			else:
				protein_seq += codon_dic[seq[i:i+3]]
		elif "N" in seq[i:i+3]:
			protein_seq += "X"

			
	return protein_seq


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
	cds_file = config_json["cdsfile"]
	pep_file = config_json["pepfile"]
	codon_file = config_json["codonfile"]
	flags_file = config_json["mutmapper"]["outdir"] + "/cds2pep-sanityflags.txt"



	codon_dic = {}
	with open(codon_file, "r") as tsvfile:
		for row in tsvfile:
			row = row.strip().split('\t')
			codon_dic[row[0]] = row[1]

	cds_seq_hash = {}
        i = 0
	for record in Seq_iO.parse(cds_file, "fasta"):
                cds_seq_hash[record.id] = str(record.seq.upper())
        	i += 1
	print "loaded %s cds sequences" % (i)


	pep_seq_hash = {}
	i = 0
	for record in Seq_iO.parse(pep_file, "fasta"):
		desc = record.description
		trs_id = desc.split("transcript:")[1].split()[0]
		pep_seq_hash[trs_id] = str(record.seq.upper())
		i += 1
	print "loaded %s pep sequences" % (i)
	

	FW = open(flags_file, "w")
	for trs_id in cds_seq_hash:
		flag = 0
		translated_pep = my_translator(cds_seq_hash[trs_id], codon_dic)
		pep_seq_hash[trs_id] = translated_pep.replace("U", "*")
		if translated_pep == pep_seq_hash[trs_id]:
			flag = 1
		FW.write("%s|%s\n" % (flag, trs_id ))
	FW.close()


if __name__ == '__main__':
        main()








