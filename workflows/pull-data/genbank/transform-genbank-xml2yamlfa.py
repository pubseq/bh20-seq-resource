#!/usr/bin/env python3
#
# Create a single YAML/FASTA for each genbank entry in GenBank XML file
#
#   transform-genbank-xml2yamlfa --out ~/tmp/pubseq file(s)
#
# Also writes a validation file in the outdir named state.json
# ----------------------------------------------------------------------

# See also directory .guix-run and README.md

import argparse
import gzip
import json
import os
import sys
import types
import xml.etree.ElementTree as ET
from utils import chunks
import genbank

parser = argparse.ArgumentParser()
parser.add_argument('--out', type=str, help='Directory to write to',
required=True)
parser.add_argument('files', nargs='+', help='file(s)')
args = parser.parse_args()

dir = args.out
if not os.path.exists(dir):
    raise Exception(f"Directory {dir} does not exist")

states = {}

for xmlfn in args.files:
    print(f"--- Reading {xmlfn}")
    try:
        with gzip.open(xmlfn, 'r') as f:
            xml = f.read().decode()
    except Exception:
        with open(xmlfn, 'r') as f:
            xml = f.read()
    tree = ET.fromstring(xml)
    for gb in tree.findall('./GBSeq'):
        valid = None
        error = None
        meta = {}
        id = gb.find("GBSeq_locus").text
        basename = dir+"/"+id
        print(f"    parsing {id}")
        try:
            valid,meta = genbank.get_metadata(id,gb)
            if valid:
                # --- write JSON
                jsonfn = basename + ".json"
                with open(jsonfn, 'w') as outfile:
                    print(f"    writing {jsonfn}")
                    json.dump(meta, outfile, indent=4)
                # --- write FASTA
                fa = basename+".fa"
                seq = genbank.get_sequence(id,gb)
                print(f"    writing {fa}")
                with open(fa,"w") as f2:
                    f2.write(f"> {id}\n")
                    f2.write(seq)
                # print(seq)
        except genbank.GBError as e:
            error = f"{e} for {id}"
            print(error,file=sys.stderr)
            valid = False
        state = {}
        state['valid'] = valid
        if error:
            state['error'] = error
        if meta['warnings']:
            state['warnings'] = meta['warnings']
        states[id] = state

statefn = dir + '/state.json'
with open(statefn, 'w') as outfile:
    print(f"    Writing {statefn}")
    json.dump(states, outfile, indent=4)
