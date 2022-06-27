import os,sys
import string
from optparse import Option_parser
import glob
import json
from Bio import Seq_iO
from Bio.Seq import Seq

__version__="4.0"
__status__ = "Dev"

###################################################
def main():
	
	usage = "\n%prog  [options]"
        parser = Option_parser(usage,version="%prog " + __version__)
	parser.add_option("-i", "--configfile", action = "store", dest = "configfile", help = "Config json file path")
	parser.add_option("-s","--srctype",action="store",dest="srctype",help="Source type")
        
	(options,args) = parser.parse_args()
        for file in ([options.configfile, options.srctype]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

	
	src_type = options.srctype
	json_file = options.configfile
	config_json = json.loads(open(json_file).read())
	genome_file = config_json["genomefile"]
        	
	output_dir = config_json["datareporter"][src_type]["outputdir"]
	input_dir = config_json["datareporter"][src_type]["inputdir"]


        project_list_file = config_json["datareporter"]["tcga"]["downloaddir"] + "/PROJECTLIST.txt"

	file_list = []
	if src_type.lower() == "tcga":
		with open(project_list_file, "r") as FR:
                	for line in FR:
                        	prj_name = line.strip()
                        	if prj_name != "" and prj_name[0:1] != "#":
					file_list += glob.glob(input_dir + prj_name + "/*/*.vcf")
	else:
		pattern = input_dir + "*.vcf"
		file_list = glob.glob(pattern)

	#Load genomic sequences
        genome_seq = {}
        for record in Seq_iO.parse(genome_file, "fasta"):
                chr_id = record.id.strip().replace("chr", "")
		genome_seq[chr_id] = record.seq.upper()

	#Initialize progress file	
	progress_file = output_dir + "/check-vcf-progress.txt"
        with open(progress_file, "w") as FW:
        	FW.write("Started writing to progress file\n")

	#Initialize bad_rows file
	bad_rows_file = output_dir + "badrows.txt"
	with open(bad_rows_file, "w") as FW:
		FW.write("")


	seen = {"chrid":{}}	

	file_count = 0
	for vcf_file in file_list:
		file_count += 1
		FR = open(vcf_file, "r")
		chr_range = map(str,range(1,23)) + ["X","Y"]
		line_count = 0
		FA1 = open(bad_rows_file, "a")
		for line in FR:
			line_count += 1
			if line.startswith('#'):
				continue
	                values = line.strip().split("\t")
				
			if src_type.lower() in ["civic"]:
				mut_id = "."
				chr_id,pos,ref,alt = values[0:4]
			else:
				chr_id,pos,mut_id,ref,alt = values[0:5]
 
			error_list = []
			if src_type.lower() in ["icgc","cosmic","clinvar"] and len(values) != 8:
				error_list.append("type1")
			if src_type.lower() == "tcga" and len(values) != 11:
				error_list.append("type2")
			if chr_id not in chr_range and not chr_id.startswith("<GL"):
				error_list.append("type3")
			if not pos.isdigit(): #Accept only integer
				error_list.append("type4")
			if not mut_id[2::].isdigit() and not mut_id[4::].isdigit() and mut_id != '.' and len(mut_id.split("t")) != 2:
				error_list.append("type5")
			if len(ref) > 1:
				error_list.append("type6")
			if len(alt) > 1:
                                error_list.append("type7")
			if ref.upper() not in ["T","A","C","G"]:
				error_list.append("type8")
			if alt.upper() not in ["T","C","G","A",",","N","-","."]:
				error_list.append("type9")
			if chr_id in genome_seq:
				if ref.upper() != genome_seq[chr_id][int(pos)-1].upper():
					error_list.append("type10")

			if len(error_list) > 0 :
				FA1.write("%s,%s,%s\n" % ("|".join(error_list), line_count, vcf_file)) 
		FA1.close()

		FA2 = open(progress_file, "a")
		FA2.write("file_index=%s  line_count=%s\n" % (file_count, line_count))
		FA2.close()
        	FR.close()


if __name__ == '__main__':
        main()





