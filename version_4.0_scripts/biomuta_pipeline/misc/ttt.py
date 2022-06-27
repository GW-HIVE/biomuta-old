in_file = "/data/projects/biomuta/downloads/v-4.0/tcga/PROJECTLIST.txt"
seen = {"one":{}, "two":{}}
with open(in_file, "r") as FR:
    for line in FR:
        c = line.strip()
        seen["one"][c] = True

in_file = "tmp/junk.csv"
with open(in_file, "r") as FR:
    for line in FR:
        parts = line.strip().replace("\"", "").split(",")
        seen["two"][parts[0]] = parts



in_file = "/data/projects/biomuta/generated/v-4.0/datareporter/doid/doid-mapping.reduced.tsv"

with open(in_file, "r") as FR:
    for line in FR:
        parts = line.split("\t")
        c = parts[2].strip()
        if c in seen["one"]:
            for do_id in parts[1].split("|"):
                if do_id in seen["two"]:
                    newrow = [c] + seen["two"][do_id]
                    print "\"%s\"" % ("\",\"".join(newrow))


