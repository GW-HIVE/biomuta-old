import csv
import json
import sys
import re


from optparse import Option_parser

__version__ = "4.0"
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





############################################################3333
def main():
	
	usage = "\n%prog [options]"
	parser = Option_parser(usage, version = "%prog" + __version__)
	parser.add_option("-i", "--configfile", action = "store", dest = "configfile" , help = "Config Json File Path")
	
	(options, args) = parser.parse_args()
	for file in ([options.configfile]):
		if not(file):
			parser.print_help()
			sys.exit(0) 

	config_file = options.configfile
	config_json = json.loads(open(config_file).read())
	
	in_file = config_json["datareporter"]["uniprot"]["input_dir"] + "/uniprot_sprot.dat"
	with open(in_file, 'r') as FR:
                row_count = 0
		ac = ""
		var_list = []
		flag = False
		for line in FR:
			row_count +=1
			prefix = line[0:2]
			if prefix == "AC":
				ac = line[3:].strip().split(";")[0]
			if line[0:17] == "OS   Homo sapiens":
				flag = True	
			if line[0:12] == "FT   VARIANT":
				parts = re.sub(r'(\s+)', "|", line).split("|")
				if parts[2] == parts[3] and parts[5] == "->":
					comment = " ".join(parts[7:])
					val = "%s|%s|%s|%s" % (parts[2], parts[4], parts[6], comment) 
					var_list.append(val)
			if line[0:2] == "//" and flag == True and len(var_list) > 0:
				var_list = sorted(set(var_list))
				for var in var_list:
					print "%s|%s" % (ac, var)
				ac = ""
				var_list = []
				flag = False

if __name__ == "__main__":
	main()



