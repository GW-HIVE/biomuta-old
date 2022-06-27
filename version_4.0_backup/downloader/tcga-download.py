import os
#import downloader
import sys
import commands
from optparse import OptionParser

__version__="1.0"
__status__ = "Dev"


###########################
def main():

    wrk_dir = "/data/projects/biomuta/downloads/v-5.0/tcga/"
    t_file = wrk_dir + "gdc-user-token.2019-07-23T17_13_24.556Z.txt"
    in_file = wrk_dir + "PROJECTLIST.txt"
    
    with open(in_file, "r") as FR:
        for line in FR:
            cancer_type = line.strip()
            m_file = wrk_dir + "manifest-files-grch38/%s.txt" % (cancer_type)
            out_dir = wrk_dir + "%s/" % (cancer_type)
            if os.path.isdir(out_dir) == False:
                cmd = "mkdir %s" % (out_dir)
                x = commands.getoutput(cmd)
            else:
                cmd = "rm -rf %s/*" % (out_dir)
                x = commands.getoutput(cmd)

            log_file = "%s.log" % (cancer_type)
            cmd = "/usr/bin/gdc-client download -t %s -m %s -d %s > %s " % (t_file, m_file, out_dir, log_file)
            x = commands.getoutput(cmd)
            





if __name__ == '__main__':
	main()
