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
       

	in_file = config_json["loader"]["outdir"] + "/do2uberon.csv"
	#Insert placeholder first
	seen = {}
	with open(in_file, 'r') as FR:
                csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
                for row in csv_grid:
			do_id1,do_id2,ro_id, uberon_id = row[0], row[1], row[2], row[3]
			if do_id1 == "-":
				continue
			if uberon_id == "-":
                                continue
			do_id1 = do_id1.split("_")[1]
			do_id2 = do_id2.split("_")[1]
			ro_id = ro_id.split("_")[1]
			uberon_id = uberon_id.split("_")[1]
			
			pair = do_id1 + ":" + uberon_id
			if pair not in seen:
				sql = "INSERT INTO biomuta_do2uberon (do_id,uberon_id) VALUES ('%s','%s') " % (do_id1, uberon_id)
				seen[pair] = True
				cur.execute(sql)
				print sql

			pair = do_id2 + ":" + uberon_id
                        if pair not in seen:
                                sql = "INSERT INTO biomuta_do2uberon (do_id,uberon_id) VALUES ('%s','%s') " % (do_id2, uberon_id)
                                seen[pair] = True
                                cur.execute(sql)
                                print sql


	DBH.commit()

if __name__ == '__main__':
        main()








