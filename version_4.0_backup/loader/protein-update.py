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
	parser.add_option("-i","--configfile",action="store",dest="configfile",help="NT file")


	(options,args) = parser.parse_args()
	for file in ([options.configfile]):
		if not (file):
			parser.print_help()
			sys.exit(0)




	config_json = json.loads(open(options.configfile, "r").read())
 	DBH = My_sQLdb.connect(host = config_json["dbinfo"]["host"], 
				user = config_json["dbinfo"]["userid"], 
				passwd = config_json["dbinfo"]["password"],
				db = config_json["dbinfo"]["dbname"])
        
	cur = DBH.cursor()
      
 
	in_file = config_json["loader"]["outdir"] + "/recnames.csv"
	i = 0
	with open(in_file, 'r') as FR:
                csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
                for row in csv_grid:
			val = row[1].replace("'", "`")
			string = "UPDATE biomuta_protein set description = '%s' WHERE canonical_ac LIKE '%s' "
			sql = string % (val,row[0] + "%")
			#print sql
			cur.execute(sql)
			i += 1
			if i%1000 == 0:
				print "Loaded %s lines " % (i)

	DBH.commit()

if __name__ == '__main__':
        main()








