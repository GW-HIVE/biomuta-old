import os,sys
import string
import csv
import json
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq


__version__="1.0"
__status__ = "Dev"


#################################
def reportProgress(progressFile, msg, openNew=False):

        if openNew == True:
                with open(progressFile, "w") as FW:
                        FW.write(msg)
        else:
                with open(progressFile, "a") as FA:
                        FA.write(msg)
        return


###################
def getAttrDict(attrString):

        attHash = {}
        attrList = attrString.strip().split(";")
        for attr in attrList:
                if attr != "":
                        name, value = attr.strip().split(" ")
                        attHash[name] = value.replace("\"", "")
        return attHash


###########################
def correctHarvestedCds(harvestedSeq, frame1, frame2):

	if frame1 == 2:
		return "a", "N" + harvestedSeq
	elif frame1 == 1:
		return "b", "NN" + harvestedSeq
	elif frame2 == 1:
                return "c", harvestedSeq[:-4] + harvestedSeq[-3:]
	elif frame2 == 2:
                return "d", harvestedSeq[:-5] + harvestedSeq[-3:]
	else:
		return "x", harvestedSeq






########################################
def loadTranscriptObjects(gtfFile, chrIdList, genomeSeqHash, progressFile):


        
        
	exonObjs = {}
        with open(gtfFile, 'rb') as tsvfile:
                tsvreader = csv.reader(tsvfile, delimiter='\t', quotechar='|')
                rowCount = 0
                for row in tsvreader:
                        if len(row) > 2:
				condList = []
				condList.append(row[0] in chrIdList)
				condList.append(row[2] in ["CDS", "stop_codon"])
				if False not in condList:
					attHash = getAttrDict(row[-1])
					obj = {
						"cds_type":row[2],
						"chr_id":row[0],
						"gene_name":attHash["gene_name"],
						"transcript_id":attHash["transcript_id"],
						"frame":int(row[7]),
						"exon_number":int(attHash["exon_number"]),
						"start":int(row[3]), "stop":int(row[4]),
						"strand":row[6]
					}
					chrId = obj["chr_id"]
					trsId = obj["transcript_id"]
					i, j = int(obj["start"]) - 1, int(obj["stop"])
                                        obj["seq"] = ""
                                        if chrId in genomeSeqHash:
                                                obj["seq"] = genomeSeqHash[chrId][i:j]
                                                if obj["strand"].strip() == "-":
                                                        tmpSeq = Seq(str(obj["seq"]))
                                                        obj["seq"] =  tmpSeq.reverse_complement()
                                        obj["seq"] = str(obj["seq"])
					if trsId not in exonObjs:
                                                exonObjs[trsId] = {"exonlist":[]}
                                        exonObjs[trsId]["exonlist"].append(obj)

                                	if rowCount%25000 == 0:
						msg = "Loaded %s exons from GTF file\n" % (rowCount)
                                        	reportProgress(progressFile, msg)
                                	rowCount += 1


	return exonObjs




###############################
def main():

	usage = "\n%prog  [options]"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--configfile",action="store",dest="configfile",help="NT file")


	(options,args) = parser.parse_args()
	for file in ([options.configfile]):
		if not (file):
			parser.print_help()
			sys.exit(0)

	configJson = json.loads(open(options.configfile, "r").read())
	genomeFile = configJson["genomefile"]

	#Load genomic sequences
	genomeSeqHash = {}
        for record in SeqIO.parse(genomeFile, "fasta"):
                chrId = record.id.strip().replace("chr", "")
		if chrId in ['12']:
			seq = record.seq.upper()
			print ">12"
			print seq
			sys.exit()
	#12	protein_coding	gene	25357723	25403870







if __name__ == '__main__':
        main()








