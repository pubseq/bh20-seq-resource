# --- Normalize data
# normalize-yamlfa.py [--yaml] --in ~/tmp/pubseq/state.json file(s)
#
# Example:
#
#    python3 ./workflows/tools/normalize-yamlfa.py -s ~/tmp/yamlfa/state.json --species ncbi_host_species.csv --specimen specimen.csv --validate

import argparse
import json
import os
import sys
import types
import normalize.mapping as mapping

parser = argparse.ArgumentParser(description="""

Normalize parameters in PubSeq JSON/YAML files. All entries in
directory are parsed using the state.json file. It is possible
to select a subset of IDs.

This tool has two modes of operation. It can validate with the
--validate switch which stops at a warning and does no rewriting.
This mode is typically used in troubleshooting.

The other mode is --rewrite which rewrites the JSON files after
making a backup (.bak) of the original. This mode updates files and
won't stop - it is used for (automated) uploads.

""")

parser.add_argument('-s','--state', type=str, help='State file (JSON) as produced by transform2yamlfa', required=True)
parser.add_argument('--species', type=str, help='Species mapping file')
parser.add_argument('--specimen', type=str, help='Optional specimen mapping file')
parser.add_argument('--validate', action='store_true', help='Validation mode - stops on warning')
parser.add_argument('--rewrite', action='store_true', help='Rewrite mode - updates files')
parser.add_argument('--yaml', action='store_true', help='Input YAML instead of JSON')
parser.add_argument('id', nargs='*', help='optional id(s)')

args = parser.parse_args()

with open(args.state) as jsonf:
    data = json.load(jsonf)

dir = os.path.dirname(args.state)
do_validate = args.validate
do_rewrite = args.rewrite

ids = args.id
if not len(ids):
    ids = list(data.keys())

species = {}
if args.species:
    with open(args.species) as f:
        for line in f:
            name,uri = line.strip().split(',')
            species[name] = uri
else:
    print("WARNING: no species mapping",file=sys.stderr)
specimen = {}
if args.specimen:
    with open(args.specimen) as f:
        for line in f:
            name,uri = line.strip().split(',')
            specimen[name] = uri
else:
    print("WARNING: no specimen mapping",file=sys.stderr)

for id in ids:
    if args.yaml:
        raise Exception("YAML not yet supported")
    fn = f"{dir}/{id}.json"
    print(f"Reading {fn}",file=sys.stderr)
    with open(fn) as f:
        rec = types.SimpleNamespace(**json.load(f))
        if do_validate:
            print(rec)
        rec.host,warning = mapping.host_species(rec.host,species)
        if warning:
            print("WARNING "+warning,file=sys.stderr)
            rec.warnings.append(warning)
        rec.sample,warning = mapping.specimen_source(rec.sample,specimen)
        if warning:
            print("WARNING "+warning,file=sys.stderr)
            rec.warnings.append(warning)
        print(rec)
        if do_validate and warning:
            print("bailing out in validation mode",file=sys.stderr)
            sys.exit(2)
        if do_rewrite:
            if not os.path.exists(fn+".bak"): # make backup the first time
                os.rename(fn,fn+".bak")
            with open(fn, 'w') as outfile:
                print(f"    Writing {fn}")
                json.dump(rec.__dict__, outfile, indent=2)
        else:
            print(rec)
