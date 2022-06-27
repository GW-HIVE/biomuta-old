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
	parser.add_option("-i","--configfile",action="store",dest="configfile",help="Config file")
	parser.add_option("-v","--viewname",action="store",dest="viewname",help="Table name")


	(options,args) = parser.parse_args()
	for file in ([options.configfile, options.viewname]):
		if not (file):
			parser.print_help()
			sys.exit(0)



	config_json = json.loads(open(options.configfile, "r").read())
 	DBH = My_sQLdb.connect(host = config_json["dbinfo"]["host"], 
				user = config_json["dbinfo"]["userid"], 
				passwd = config_json["dbinfo"]["password"],
				db = config_json["dbinfo"]["dbname"])

        cur = DBH.cursor()
	sql = config_json["view_tables"][options.viewname]
	cur.execute(sql)
	DBH.commit()


if __name__ == '__main__':
        main()








