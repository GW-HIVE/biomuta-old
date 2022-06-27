import json
import requests
import csv
import sys
import yaml


"""
Pulls dbsnp.tsv file from ebi api for protein list in proteins.dat (should update to take generic user input.
Pulls in batches of 100 so as to not overwhelm the api and writes the json to batch files for later parsing if needed
These json request are then filtered and converted into dbsnp.tsv
Updated api call based on email with A.N. from ebi on 2019-May-16
Run this in a clean directory containing only this script and proteins.dat

If the connection was lost before this script finished (it can take over a day for the proteins in glygen) run again in same directory with system argument 'restart'
ex. python create_dbSNP_table_uniprot.py restart
this will read through the batches to see what was already pulled and added to dbsnp.tsv and will continue from where it left off
"""


#idmappingFile = "/data/projects/glygen/generated/datasets/reviewed/human_protein_idmapping.csv"
protListFile  = "proteins.dat"

class VAR():
    header = "\t".join(["chrom", "pos", "id", "ref", "alt", "frequency", "geneinfo", "accession", "pos_in_uniprot", "ref_aa", "alt_aa", "in_dbSNP"]) #see __repr__() and row()
    def __init__(self):
        self.accession = ""
        self.geneinfo  = ""
        self.ref_aa     = ""
        self.alt_aa     = ""
        self.frequency = ""
        self.pos_in_uniprot = ""
        self.chrom = ""
        self.pos = ""
        self.ref = ""
        self.alt = ""
        self.iD  = ""
        self.in_dbSNP = "False"
        self.has_rsID = False

    def valid(self):
        if not(self.has_rsID):
            return False
        try:
            int(self.chrom)
        except:
            return False

        return True

    def __repr__(self):
        return "\t".join([self.chrom, self.pos, self.iD, self.ref, self.alt, self.frequency, self.geneinfo, self.accession, self.pos_in_uniprot, self.ref_aa, self.alt_aa, self.in_dbSNP]) #see self.header and self.row

    def row(self):
        return [self.chrom, self.pos, self.iD, self.ref, self.alt, self.frequency, self.geneinfo, self.accession, self.pos_in_uniprot, self.ref_aa, self.alt_aa, self.in_dbSNP] #see header and self.__repr__



def get_protein_list_from_idmapping(idmappingFile):
    proteins = set(())
    with open(idmappingFile, 'r') as handle:
        next(handle) #skip header
        protCSV = csv.reader(handle, delimiter=",")
        for row in protCSV:
            proteins.add(row[0].split("-")[0]) #remove isoform number
        return proteins


def get_protein_list(protFil):
    proteins = set(())
    with open(protFil, 'r') as handle:
        next(handle) #skip header
        for row in handle:
            proteins.add(row.rstrip()) #remove isoform number
        return proteins



def restart_protein_list(protSet, downloadFile):
    #reads through download file, since download  download file is writter after variations we can pick up with the next batch
    countWritten = 0
    with open(downloadFile, 'r') as rstrtHand:
        count = 0
        next(rstrtHand) #skip header
        for line in rstrtHand:
            spl = line.split("\t")
            count += 1
            protSet.remove(spl[0])
            batch = int(spl[1])
    if not((count % 100) == 0):
        print "Warning: Writing of batch {} did not finish. Will have duplicated entries in batch.json files and variations file.".format(batch)
    return protSet, batch+1 #start at batch after one that is recoreded to have finished


def run_request(reqList):
    requestURL = "https://www.ebi.ac.uk/proteins/api/variation?size=-1&accession={}".format("%2C".join(reqList))
    #requestURL = "https://www.ebi.ac.uk/proteins/api/variation?size=-1&accession=P31689%2CP10144%2CO95298"
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.text



def parse_genomic_loc(genLoc):
    #ex NC_000014.9:g.24634144A>G
    if "ins" in genLoc:
        return "","","",""
    spl = genLoc.split(":")
    chrom = str(int(spl[0].split('.')[0][-3:]))  #NC_000014.9, takes 014 converts to int to remove 0 then back to string
    posRefAltStr = spl[1].split(".")[1] #g.24634144A>G  -> 24634144A>G
    posRefStr, alt = posRefAltStr.split(">")
    pos = ''
    ref = ''
    for character in posRefStr:   #fine for 24634144A, will break for malformed geneLocations such as 123A4G>CC, will give pos =  1234, and ref  = AG
        if character.isdigit():
            pos += character
        else:
            ref += character
    return chrom, pos, ref, alt
    



def parse_json_to_dbSNP(jsonResponse):
#chr     pos     id      ref     alt     caf     sao     geneinfo        vc      nsf     nsm     nsn     syn     trsid   pos_in_trs      pep_id  pos_in_pep      uniprot_isoform_ac      pos_in_uniprot_isoform  ref_aaalt_aa
#1       69096   rs1171757348    G       A                       OR4F5:79501     SNV                             SYN     ENST00000335137 6       ENSP00000334393 2       Q8NH21-1        2       V       V
    varList = []
    for protein in jsonResponse:
        try:
            accession = protein["accession"]
        except:
            pass

        try:
            geneinfo  = protein["geneName"]
        except:
            pass

        for var in protein["features"]:
            tempVar = VAR()
            tempVar.accession = accession
            try:
                tempVar.geneinfo  = geneinfo
            except:
                pass                            
            try:
                tempVar.ref_aa    = var["wildType"]    
            except:
                pass

            try:
                tempVar.alt_aa    = var["alternativeSequence"]
            except:
                pass

            try: #not all entries have frequency
                tempVar.frequency = var["frequency"]
            except:
                pass
            
            try:
                tempVar.pos_in_uniprot = var["begin"]
            except:
                pass
            try:
                tempVar.chrom, tempVar.pos, tempVar.ref, tempVar.alt = parse_genomic_loc(var["genomicLocation"]) 
            except:
                pass
            try:
                for xRef in var["xrefs"]:   #multiple xrefs have same rsID, will take first entry containing rsID regardless of source "ESP", "Ensembl", "ExAC", etc.
                    if (xRef["id"][0:2] == "rs"):
                        tempVar.iD = xRef["id"]
                        tempVar.has_rsID = True
                    if (xRef["name"] == "dbSNP"):
                        tempVar.in_dbSNP = "True"
            except:
                pass
            varList.append(tempVar)
    return varList

def main():
    restart = False
    if ('restart' in sys.argv):
        print "Warning: When restarting, variations and batch.json output may be duplicated for last batch completed."
        restart = True
    #proteins = get_protein_list_from_idmapping(idmappingFile) 
    proteins = get_protein_list("proteins.dat")
    print "Pulling variations for {} proteins.".format(len(proteins))
    zCount = (len(str(len(proteins)-1))) - 2 #numer of digits for zero pad formatting, (-2 because we will have 100 proteins per request batch)
    reqList = []
    batch = 0

    if not(restart):
        #write header to dbSNP.tsv - this will overwrite
        with open("dbSNP.tsv", "w") as dbSNPOut:
            dbSNPOut.write(VAR.header)                       
            dbSNPOut.write("\n")
    
        #write header to dbSNP_download.log - this will overwrite
        with open("dbSNP_download.log", 'w') as log: #record which protein download is in which batch file
            log.write('"accession"\t"batch"\n')
    
    else:
        proteins, batch = restart_protein_list(proteins, "dbSNP_download.log")

    #batch proteins into sets of 100 for api call
    while (len(proteins) > 0):
        if (len(reqList) == 100): #ebi protein variations api can only accept 100 entries at a time  
            print "---------------------------------------------"
            print "Pulling Response for batch #{}".format(str(batch))
            textResponse = run_request(reqList) #get response
            print "Uploading to Dictionary"
            jsonResponse = yaml.safe_load(textResponse) #create dictionary, avoid unicode

            #write json to a batch download file just in case it is needed laer
            print "Saving JSON"
            with open(str(batch).zfill(zCount)+"_batch.json", "w") as jsonOut:
                jsonOut.write(json.dumps(jsonResponse, indent=4, sort_keys=True)) #pretty formatting

            print "Parsing Response"
            #jsonResponse = json.loads(textResponse) #create dictionary 
            varList = parse_json_to_dbSNP(jsonResponse) #parse json to get variations

            #write to dbSNP.tsv
            print "Writing dbSNP"
            with open("dbSNP.tsv", "a") as dbSNPOut:
                dbSNPCSV = csv.writer(dbSNPOut, delimiter="\t")
                for var in varList:
                    if var.valid():
                        dbSNPCSV.writerow(var.row())


            #record which protein is in which json batch
            print "Updating Download Log"
            with open("dbSNP_download.log", 'a') as log: #record which protein download is in which batch file
                log.write("\n".join([x+"\t"+ str(batch).zfill(zCount) for x in reqList]))
                log.write("\n")
            reqList = []
            batch += 1    
   
        #if not at api max, add more accession to the reqList
        else:
            reqList.append(proteins.pop())



    #capture last batch of proteins to download (would not have been pulled in while loop unless #prots %100 == 0!
    #TODO make these steps a function so we do not need to update both here and in while loop
    if (len(reqList) > 0):
        print "---------------------------------------------"
        print "Pulling Response for batch #{}".format(str(batch))
        textResponse = run_request(reqList) #get response
        print "Uploading to Dictionary"
        jsonResponse = yaml.safe_load(textResponse) #create dictionary, avoid unicode

        #write json to a batch download file just in case it is needed laer
        print "Saving JSON"
        with open(str(batch).zfill(zCount)+"_batch.json", "w") as jsonOut:
            jsonOut.write(json.dumps(jsonResponse, indent=4, sort_keys=True)) #pretty formatting

        print "Parsing Response"
        #jsonResponse = json.loads(textResponse) #create dictionary
        varList = parse_json_to_dbSNP(jsonResponse) #parse json to get variations

        #write to dbSNP.tsv
        print "Writing dbSNP"
        with open("dbSNP.tsv", "a") as dbSNPOut:
            dbSNPCSV = csv.writer(dbSNPOut, delimiter="\t")
            for var in varList:
                if var.valid():
                    dbSNPCSV.writerow(var.row())


        #record which protein is in which json batch
        print "Updating Download Log"
        with open("dbSNP_download.log", 'a') as log: #record which protein download is in which batch file
            log.write("\n".join([x+"\t"+ str(batch).zfill(zCount) for x in reqList]))
            log.write("\n")
        #reqList = []
        #batch += 1



main()
