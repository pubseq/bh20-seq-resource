#!/usr/bin/env python3
#
# - bulk download genbank data and matadata, preparing the FASTA and
#   the YAML files
#
# update-from-genbank.py --max 10 --skip ids.txt --outdir ~/tmp/genbank
#
# See directory .guix-run and README.md

BATCH_SIZE=5000

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--max', type=int, help='Max queries', required=False)
parser.add_argument('--skip', type=str, help='File with ids to skip, 1 id per line', required=False)
parser.add_argument('--outdir', type=str, help='Output directory', required=True)
args = parser.parse_args()

from Bio import Entrez
Entrez.email = 'another_email@gmail.com' # FIXME

import xml.etree.ElementTree as ET
import json
import os
import requests

from datetime import date, datetime
from dateutil.parser import parse

import sys
# sys.path.append('../')
from utils import is_integer, chunks, check_and_get_ontology_dictionaries

num_ids_for_request = 100 # batch calls
min_acceptable_collection_date = datetime(2019, 12, 1)

outdir = args.outdir

today_date = date.today().strftime("%Y.%m.%d")
path_ncbi_virus_accession = 'sequences.{}.acc'.format(today_date)

if not os.path.exists(outdir):
    raise Exception(f"Output directory {outdir} does not exist!")

skip = set()
if args.skip:
    with open(args.skip) as f:
        content = f.readlines()
        for line in content:
            skip.add(line.strip())

print(f"Skip size is {len(skip)}",file=sys.stderr)

# ----------------------------------------------------------------------
"""
Download section for genbank XML
"""

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
            Entrez.esearch(db='nuccore', term=term, idtype='acc',
                           retstart=retstart, retmax=BATCH_SIZE)
        )
        idlist = record['IdList']
        new_ids = set(idlist)
        num_read = len(new_ids)
        print(num_read,":",idlist[0],file=sys.stderr)
        retstart += num_read
        new_ids.difference_update(skip) # remove skip ids
        new_ids = set([id for id in new_ids if id[:2] not in PREFIX])
        ids.update(new_ids)             # add to total set
        print(f"Term: {term} --> #{len(new_ids)} new IDs ---> Total unique IDs #{len(ids)})",file=sys.stderr)
        if args.max and len(ids) > args.max:
            print(f"Stopping past #{args.max} items",file=sys.stderr)
            break

for id in ids:
    print(id)

sys.exit(2)

with open(path_ncbi_virus_accession) as f:
    tmp_list = [line.strip('\n') for line in f]

new_ids_set = set(tmp_list)
if len(accession_to_consider_set) > 0:
    new_ids_set = new_ids_set.intersection(accession_to_consider_set)

new_ids = len(new_ids_set.difference(id_set))
id_set.update(new_ids_set)

print('DB: NCBI Virus', today_date, '-->', new_ids, 'new IDs from', len(tmp_list), '---> Total unique IDs:', len(id_set))

id_set = id_set.difference(accession_to_ignore_set)
print('There are {} missing IDs to download.'.format(len(id_set)))

os.makedirs(outdir)
for i, id_x_list in enumerate(chunks(list(id_set), num_ids_for_request)):
    path_metadata_xxx_xml = os.path.join(outdir, 'metadata_{}.xml'.format(i))
    print('Requesting {} ids --> {}'.format(len(id_x_list), path_metadata_xxx_xml))

    with open(path_metadata_xxx_xml, 'w') as fw:
        fw.write(
            Entrez.efetch(db='nuccore', id=id_x_list, retmode='xml').read()
        )
