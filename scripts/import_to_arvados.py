import os
import subprocess
import glob
import sys

os.chdir(os.environ["TMPDIR"])
os.symlink(sys.argv[2], "dict_ontology_standardization")
subprocess.run(sys.argv[1])

os.chdir("fasta_and_yaml")
fasta_files = glob.glob("*.fasta")

for f in fasta_files:
    subprocess.run(["bh20-seq-uploader", f, "%s.yaml" %f[:-6]])
