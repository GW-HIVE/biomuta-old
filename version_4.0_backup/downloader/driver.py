import downloader
import sys
from optparse import OptionParser

__version__="1.0"
__status__ = "Dev"


###########################
def main():
	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--configfile",action="store",dest="configfile",help="filepath for config file")
	parser.add_option("-s","--source",action="store",dest="source",help="download source")
	(options,args) = parser.parse_args()
	
	for file in ([options.configfile, options.source]):
	        if not (file):
        	        parser.print_help()
	                sys.exit(0)
	
	source = options.source
	json_file = options.configfile
	
	
	# Driver
	rs = downloader.DataResource(json_file)
	rs.download(source)
	rs.checkmd5(source)

if __name__ == '__main__':
	main()
