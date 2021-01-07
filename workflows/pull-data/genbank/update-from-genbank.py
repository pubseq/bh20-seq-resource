#!/usr/bin/env python3
#
# bulk download genbank data and matadata, preparing the FASTA and the
# YAML files
#
#   update-from-genbank.py --max 10 --ids ids.txt --out ~/tmp/genbank-xml
#
# See also directory .guix-run and README.md

import argparse
import gzip
import os
import sys
from utils import chunks

from Bio import Entrez

Entrez.email = 'another_email@gmail.com'  # FIXME

BATCH = 100

parser = argparse.ArgumentParser()
parser.add_argument('--ids', type=str, help='File with ids to fetch, 1 id per line', required=True)
parser.add_argument('--out', type=str, help='Directory to write to', required=True)
parser.add_argument('--max', type=int, help='Max queries', required=False)

args = parser.parse_args()

ids = set()
with open(args.ids) as f:
    ids.update([line.strip() for line in f])

dir = args.out
if not os.path.exists(dir):
    raise Exception(f"Directory {dir} does not exist")

request_num = BATCH
if args.max:
    request_num = min(BATCH, args.max)

for num_chunk, ids_chunk in enumerate(chunks(list(ids), request_num)):
    xmlfn = os.path.join(dir, f"metadata_{num_chunk}.xml.gz")
    print(f"Fetching {xmlfn} ({num_chunk * request_num})", file=sys.stderr)

    with gzip.open(xmlfn, 'w') as f:
        f.write(Entrez.efetch(db='nuccore', id=ids_chunk, retmode='xml').read().encode())

    if args.max and num_chunk * request_num >= args.max:
        break
