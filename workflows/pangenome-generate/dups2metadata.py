import sys
import re

md = open(sys.argv[1], "rt")
for d in md:
    sys.stdout.write(d)

if len(sys.argv) < 3:
    exit(0)

sameseqs = open(sys.argv[2], "rt")
for d in sameseqs:
    g = re.match(r"\d+\t(.*)", d)
    sp = g.group(1).split(",")
    for n in sp[1:]:
        sys.stdout.write("<%s> <http://biohackathon.org/bh20-seq-schema/has_duplicate_sequence> <%s> .\n" % (n.strip(), sp[0].strip()))
