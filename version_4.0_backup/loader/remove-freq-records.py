import os,sys
import string
import csv
import json
from optparse import Option_parser
import My_sQLdb



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

	config_json = json.loads(open(options.configfile, "r").read())
 	DBH = My_sQLdb.connect(host = config_json["dbinfo"]["host"], 
				user = config_json["dbinfo"]["userid"], 
				passwd = config_json["dbinfo"]["password"],
				db = config_json["dbinfo"]["dbname"])

	progress_file = config_json["loader"]["outdir"] + "/deleterecords-freq-progress.txt"
        cur = DBH.cursor()

	report_progress(progress_file, "Started deleting records", True)
	for i in xrange(1,500):
		sql = "DELETE FROM biomuta_mutation_freq WHERE data_src = '%s' LIMIT 100000 " % (src_type)
		cur.execute(sql)
		DBH.commit()
		report_progress(progress_file, "%s batch deleted " % (i))


if __name__ == '__main__':
        main()








