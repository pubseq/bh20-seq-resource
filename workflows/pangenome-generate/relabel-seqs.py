import os
import json

reads = []
b = 1
while os.path.exists("block%i" % b):
    with open("block%i" % b) as f:
        reads.extend(json.load(f))
    b += 1

with open("subjects") as f:
    subjects = json.load(f)

relabeled_fasta = open("relabeledSeqs.fasta", "wt")
original_labels = open("originalLabels.ttl", "wt")

for i, r in enumerate(reads):
    with open(r["path"], "rt") as fa:
        label = fa.readline()
        original_labels.write("<%s> <http://biohackathon.org/bh20-seq-schema/original_fasta_label> \"%s\" .\n" % (subjects[i], label[1:].strip().replace('"', '\\"')))
        relabeled_fasta.write(">"+subjects[i]+"\n")
        data = fa.read(8096)
        while data:
            relabeled_fasta.write(data)
            endswithnewline = data.endswith("\n")
            data = fa.read(8096)
        if not endswithnewline:
            relabeled_fasta.write("\n")
