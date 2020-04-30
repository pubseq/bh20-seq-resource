import os
import subprocess
import glob
import sys

os.chdir(os.environ["TMPDIR"])
subprocess.run(sys.argv[1])

os.chdir("fasta_and_yaml")
fasta_files = glob.glob("*.fasta")

for f in fasta_files:
    subprocess.run(["bh20-seq-uploader", f, "%s.yaml" %f[:-6]])
