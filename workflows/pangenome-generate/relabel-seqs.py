import os
import json
import sys

def readitems(stem):
    items = []
    b = 1
    while os.path.exists("%s%i" % (stem, b)):
        with open("%s%i" % (stem, b)) as f:
            items.extend(json.load(f))
        b += 1
    return items

reads = readitems("block")
subjects = readitems("subs")

relabeled_fasta = open("relabeledSeqs.fasta", "wt")
original_labels = open("originalLabels.ttl", "wt")

blacklist = set()
if len(sys.argv) > 1:
    with open(sys.argv[1]) as bl:
        for l in bl:
            blacklist.add(l.strip())

for i, r in enumerate(reads):
    with open(r["path"], "rt") as fa:
        label = fa.readline().strip()
        original_labels.write("<%s> <http://biohackathon.org/bh20-seq-schema/original_fasta_label> \"%s\" .\n" % (subjects[i], label[1:].replace('"', '\\"')))
        skip = (subjects[i] in blacklist or label[1:] in blacklist)
        if skip:
            original_labels.write("<%s> <http://biohackathon.org/bh20-seq-schema/excluded_from_graph> \"true\"^^<http://www.w3.org/2001/XMLSchema#boolean> .\n" % (subjects[i]))
        if not skip:
            relabeled_fasta.write(">"+subjects[i]+"\n")
        data = fa.read(8096)
        while data:
            if not skip:
                relabeled_fasta.write(data)
            endswithnewline = data.endswith("\n")
            data = fa.read(8096)
        if not skip and not endswithnewline:
            relabeled_fasta.write("\n")
