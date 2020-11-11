from rdflib import Graph
import sys
import subprocess
g = Graph()
g.parse(sys.argv[1], format="nt")
res = g.query(sys.argv[3])
for r in res:
    subprocess.run(["samtools", "faidx", sys.argv[2], r[0]])
