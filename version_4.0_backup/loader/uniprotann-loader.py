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


	(options,args) = parser.parse_args()
	for file in ([options.configfile]):
		if not (file):
			parser.print_help()
			sys.exit(0)




	ann_type = "uniprot"
	config_json = json.loads(open(options.configfile, "r").read())
 	DBH = My_sQLdb.connect(host = config_json["dbinfo"]["host"], 
				user = config_json["dbinfo"]["userid"], 
				passwd = config_json["dbinfo"]["password"],
				db = config_json["dbinfo"]["dbname"])
        
	cur = DBH.cursor()

       
	in_file = config_json["loader"]["outdir"] + "/uniprot-ann.csv"
	with open(in_file, 'r') as FR:
                csv_grid = csv.reader(FR, delimiter=',', quotechar='"')
                for row in csv_grid:
			canonical_ac, ann_name, start_pos, end_pos, ann_value = row[0],row[1],row[2],row[3],row[4]
			ann_value = ann_value.replace("'", "`")
			if "-" in [row[2],row[3]]:
				continue		
	
			sql = "INSERT INTO biomuta_protein_ann "
			sql += "(canonical_ac, start_pos, end_pos, ref, alt, ann_type, ann_name,ann_value) "
			sql += "VALUES ('%s',%s,%s,'%s','%s','%s', '%s','%s')" % (canonical_ac, start_pos, end_pos,'', '', ann_type,ann_name,ann_value)
			
			print sql
			cur.execute(sql)

	
	DBH.commit()

if __name__ == '__main__':
        main()








