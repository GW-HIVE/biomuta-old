import os,sys
import string
import commands
from optparse import OptionParser
import glob
import json
import csv


__status__ = "Dev"
__version__ = "4.0"








##############################################
def main():

    usage = "\n%prog  [options]"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i", "--configfile", action = "store", dest = "configfile", help = "Config file")
 

    (options,args) = parser.parse_args()
    for file in ([options.configfile]):
        if not (file):
            parser.print_help()
            sys.exit(0)

    json_file = options.configfile
    config_json = json.loads(open(json_file).read())

    chr_list = map(str,range(1,23)) + ["X","Y"]
    #chr_list = ["1"]
 

    #Load TCGA to DO mapping
    map_file = config_json["datareporter"]["doid"]["outputdir"] + "/doid-mapping.reduced.reduced.csv"
    tcga2do = {}
    lines = open(map_file, "r").read().split("\n")
    for line in lines[1:]:
        if line.strip() != "":
            parts = line.strip().replace("\"", "").split(",")
            tcga2do[parts[0]] = {"doid":parts[1], "doname":parts[2]}


    out_file = "outdir/mutid-map.final.csv"
    FW = open(out_file, "w")
    
    row = ["chr_id","chr_pos","ref_nt","alt_nt","mut_xref_id", "uniprot_canonical_ac","aa_pos","ref_aa","alt_aa"]
    row += ["cancer_type","patients_positive", "patients_tested","mut_freq", "source", "do_id", "do_name"]
    FW.write("%s\n" % (",".join(row)))

    for chr_id_out in chr_list:
        count_dict = {}
        data_frame = {}
        src_list = []
        file_list = glob.glob("outdir/mutid-icgc-chr%s.csv" % (chr_id_out))
        file_list += glob.glob("outdir/mutid-tcga-chr%s.2.csv" % (chr_id_out))
        for in_file in file_list:
            file_name = in_file.split("/")[-1]
            src_type = file_name.split("-")[1]
            src_list.append(src_type)
            with open(in_file, 'r') as FR:
                row_count = 0
                for line in FR:
                    row_count += 1
                    row = line.split(",")
                    mut_freq = "%5.4f" % (float(row[10])/float(row[11]))
                    #print row
                    cancer_type = row[9]
                    if cancer_type in tcga2do:
                        do_id = tcga2do[cancer_type]["doid"]
                        do_name = tcga2do[cancer_type]["doname"]
                        FW.write("%s,%s,%s,%s,%s\n" % (line.strip(), mut_freq, src_type, do_id, do_name))
    FW.close()



if __name__ == '__main__':
	main()


