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




###############################
def main():

	usage = "\n%prog  [options]"
        parser = Option_parser(usage,version="%prog " + __version__)
        parser.add_option("-i", "--configfile", action = "store", dest = "configfile", help = "Config json file path")
        parser.add_option("-s","--srctype",action="store",dest="srctype",help="Source type")

        (options,args) = parser.parse_args()
        for file in ([options.srctype]):
                if not (file):
                        parser.print_help()
                        sys.exit(0)

        src_type = options.srctype
        json_file = options.configfile
        config_json = json.loads(open(json_file).read())



	if src_type != "tcga":
		print "This step applies for TCGA data_src only!"
		sys.exit()
	
	chr_list = map(str,range(1,23)) + ["X","Y"]
	
	progress_file = config_json["datareporter"][src_type]["output_dir"] + "/partition-mutlist-progress.txt"
	mutlist_file = config_json["datareporter"][src_type]["output_dir"]
	mutlist_file += "/mutlist.chrall.csv"
	
	report_progress(progress_file, "Started partitioning", True)	
	for chr_id in chr_list:
		out_file = config_json["datareporter"][src_type]["output_dir"] + "/mutlist.chr"+chr_id+".csv"
		FW = open(out_file, "w")
		with open(mutlist_file, 'rb') as tsvfile:
			tsvreader = csv.reader(tsvfile, delimiter=',', quotechar='"')
			row_count = 0
			for row in tsvreader:
				row_count += 1
				if row_count == 1:
                                     	FW.write("%s\n" % (",".join(row)))
					continue
				if row[0] != chr_id:
					continue 
				FW.write("%s\n" % (",".join(row)))
		report_progress(progress_file, "finished %s" % (chr_id))
		FW.close()




if __name__ == '__main__':
        main()








