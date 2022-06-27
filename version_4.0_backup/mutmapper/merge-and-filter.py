import os,sys
import string
import csv
import json
from optparse import Option_parser
from Bio import Seq_iO
from Bio.Seq import Seq


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



##########################################
def mask_regions(gtf_file, seq_hash, chr_id):

        with open(gtf_file, 'r') as FR:
                reader = csv.reader(FR, delimiter='\t', quotechar='|')
                row_count = 0
                n = 0
                for row in reader:
                        row_count += 1
                        if row[0][0:1] == "#" or len(row) < 2:
                                continue
                        if row[0] != chr_id:
                                continue
                        feat_type,start_pos, end_pos = row[2], int(row[3]), int(row[4])
                        if feat_type in ["CDS", "stop_codon"]:
                                for seq_index in xrange(start_pos-1, end_pos):
                                        seq_hash[chr_id][seq_index] = "1"
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
	gtf_file = config_json["gtffile"]
	genome_file = config_json["genomefile"]
	merge_file = config_json["mutmapper"]["outdir"] + "/mutlist-merged.csv"
	progress_file = config_json["mutmapper"]["outdir"] + "/merge-progress.txt"
	stat_file = config_json["mutmapper"]["outdir"] + "/merge-stat.csv"

        mut_index = 1
	chr_list = map(str,range(1,23)) + ["X","Y"]


	field_list = ["index", "chr","pos","ref","alt","codingflag"]
	FW = open(merge_file, "w")
	FW.write("%s\n" % (",".join(field_list)))

	report_progress(progress_file, "Progress report", True)
	report_progress(stat_file, "chr_id,all,passed,unique,coding,noncoding", True)
	t1, t2, t3, t4, t5 = 0, 0, 0, 0, 0
	for chr_id in chr_list:
        	chr_file = config_json["chrdir"] + "/chr"+chr_id+".fa"
		seq_hash = {}
		report_progress(progress_file, "Started loading sequence for chr %s" % (chr_id))
        	for record in Seq_iO.parse(chr_file, "fasta"):
                	seq_hash[chr_id] = list(str(record.seq.upper()))
		report_progress(progress_file, "done.")

        	report_progress(progress_file, "Started masking regions for chr %s" % (chr_id))
        	mask_regions(gtf_file, seq_hash, chr_id)
        	report_progress(progress_file, "done!")

		seen = {"mut":{}}
		n1, n2, n3, n4, n5 = 0, 0, 0, 0, 0
        	for src_name in config_json["vcfsources"]:
			mutlist_file = config_json["datareporter"][src_name]["output_dir"]
			mutlist_file += "/mutlist.chr"+chr_id+".csv"
			if os.path.exists(mutlist_file) == False:
				continue
			report_progress(progress_file, "Now parsing %s" % (mutlist_file))
			with open(mutlist_file, 'rb') as tsvfile:
                        	tsvreader = csv.reader(tsvfile, delimiter=',', quotechar='"')
                        	row_count = 0
				for row in tsvreader:
                                	row_count += 1
					if row_count == 1:
                                        	continue
                                	if row[0] != chr_id:
                                                continue
					n1 += 1
					if row[4] == "False" or row[3] not in ["A", "C", "G", "T"]:
                                        	continue
					n2 += 1
					combo_id = ",".join([row[0], row[1], row[2], row[3]])
        				if n3%100000 == 0:
						msg = "%s,%s,%s mutations" % (chr_id, src_name, n3)
                                                report_progress(progress_file, msg)
					if combo_id not in seen["mut"]:
                                        	n3 += 1
                                        	seen["mut"][combo_id] = True
                                        	seq_index = int(row[1]) - 1
                                        	coding_flag = "noncoding"
						if seq_hash[chr_id][seq_index] == "1":
                                               		coding_flag = "coding"
							n4 += 1
						else:
							n5 += 1 
						FW.write("%s,%s,%s\n" % (mut_index,combo_id,coding_flag))
						mut_index += 1
		report_progress(stat_file, "%s,%s,%s,%s,%s,%s" % (chr_id,n1,n2,n3,n4,n5))
		t1, t2, t3, t4, t5 = t1 + n1, t2 + n2, t3 + n3, t4 + n4, t5 + n5
	report_progress(stat_file, "%s,%s,%s,%s,%s,%s" % ("total",t1,t2,t3,t4,t5))
	FW.close()
	

 



if __name__ == '__main__':
        main()








