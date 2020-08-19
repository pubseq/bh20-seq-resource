import sys

md = open(sys.argv[1], "rt")
for d in md:
    print(d)

if len(sys.argv) < 3:
    exit(0)

sameseqs = open(sys.argv[2], "rt")
for d in sameseqs:
    logging.warn(d)
    g = re.match(r"\d+\t(.*)", d)
    logging.warn("%s", g.group(1))
    sp = g.group(1).split(",")
    for n in sp[1:]:
        print("<%s> <http://biohackathon.org/bh20-seq-schema/has_duplicate_sequence> <%s> ." % (n.strip(), sp[0].strip()))
