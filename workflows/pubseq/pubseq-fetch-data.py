#!/usr/bin/env python3

import argparse
import json
import os
import requests
import sys
import time

parser = argparse.ArgumentParser(description="""

Fetch metadata (JSON) from PubSeq and optionally the FASTA files.  IDs
can be passed in on the command line or in a file.

""")
parser.add_argument('--fasta', action='store_true', help='Also fetch FASTA records')
parser.add_argument('--out', type=str, help='Directory to write to',
required=True)
parser.add_argument('--ids', type=str, help='File with ids', required=False)
parser.add_argument('id', nargs='*', help='id(s)')
args = parser.parse_args()

dir = args.out
if not os.path.exists(dir):
    raise Exception(f"Directory {dir} does not exist")

ids = args.id
if (len(ids)==0):
    print(f"Reading {args.ids}")
    with open(args.ids) as f:
        ids = [ l.strip() for l in f.readlines() ]

for id in ids:
    print(id)
    jsonfn = dir+"/"+id+".json"
    if not os.path.exists(jsonfn):
        count = 0
        r = requests.get(f"http://covid19.genenetwork.org/api/sample/{id}.json")
        while not r:
            count += 1
            if count>10: raise Exception(f"Can not find record for {id}")
            time.sleep(15)
            r = requests.get(f"http://covid19.genenetwork.org/api/sample/{id}.json")
        m_url = r.json()[0]['metadata']
        mr = requests.get(m_url)
        with open(dir+"/"+id+".json","w") as outf:
            outf.write(mr.text)
        if args.fasta:
            fastafn = dir+"/"+id+".fa"
            if os.path.exists(fastafn): continue
            fa_url = r.json()[0]['fasta']
            fr = requests.get(fa_url)
            with open(fastafn,"w") as outf:
                outf.write(fr.text)

