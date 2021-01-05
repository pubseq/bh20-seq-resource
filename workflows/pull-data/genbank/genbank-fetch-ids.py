#!/usr/bin/env python3
#
# Find all genbank IDs
#
#   genbank-fetch-ids.py --max 1000 --skip pubseq_ids.txt
#
# See also directory .guix-run and README.md

import argparse
import sys

from Bio import Entrez

parser = argparse.ArgumentParser()
parser.add_argument('--max', type=int, help='Max queries', required=False)
parser.add_argument('--skip', type=str, help='File with ids to skip, 1 id per line', required=False)
args = parser.parse_args()

BATCH_SIZE = 5000

Entrez.email = 'another_email@gmail.com'  # FIXME

skip = set()
if args.skip:
    with open(args.skip) as f:
        content = f.readlines()
        for line in content:
            skip.add(line.strip())

print(f"Skip size is {len(skip)}", file=sys.stderr)

# Try to search several strings
TERMS = ['SARS-CoV-2', 'SARS-CoV2', 'SARS CoV2', 'SARSCoV2', 'txid2697049[Organism]']

# Remove mRNAs, ncRNAs, Proteins, and predicted models (more information here: https://en.wikipedia.org/wiki/RefSeq) starting with
PREFIX = ['NM', 'NR', 'NP', 'XM', 'XR', 'XP', 'WP']

ids = set()
for term in TERMS:
    num_read = BATCH_SIZE
    retstart = 0

    while num_read == BATCH_SIZE:
        record = Entrez.read(
            Entrez.esearch(db='nuccore', term=term, idtype='acc', retstart=retstart, retmax=BATCH_SIZE)
        )

        idlist = record['IdList']
        new_ids = set(idlist)
        num_read = len(new_ids)
        retstart += num_read

        print(num_read, ":", idlist[0], file=sys.stderr)

        new_ids.difference_update(skip)  # remove skip ids
        new_ids = set([id for id in new_ids if id[:2] not in PREFIX])
        ids.update(new_ids)  # add to total set

        print(f"Term: {term} --> #{len(new_ids)} new IDs ---> Total unique IDs #{len(ids)}", file=sys.stderr)

        if args.max and len(ids) > args.max:
            print(f"Stopping past #{args.max} items", file=sys.stderr)
            break

for id in ids:
    print(id)
