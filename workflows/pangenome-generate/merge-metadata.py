import re
import schema_salad.schema
import schema_salad.jsonld_context
import json
import sys
import os
import logging

metadataSchema = sys.argv[1]
originalLabels = sys.argv[2]
dups = None
if len(sys.argv) == 4:
    dups = sys.argv[3]

def readitems(stem):
    items = []
    b = 1
    while os.path.exists("%s%i" % (stem, b)):
        with open("%s%i" % (stem, b)) as f:
            items.extend(json.load(f))
        b += 1
    return items

metadata = readitems("block")
subjects = readitems("subs")

(document_loader,
 avsc_names,
 schema_metadata,
 metaschema_loader) = schema_salad.schema.load_schema(metadataSchema)

for i, m in enumerate(metadata):
    doc, metadata = schema_salad.schema.load_and_validate(document_loader, avsc_names, m["path"], False, False)
    doc["id"] = subjects[i]
    g = schema_salad.jsonld_context.makerdf(subjects[i], doc, document_loader.ctx)
    print(g.serialize(format="ntriples").decode("utf-8"))

if dups:
    sameseqs = open(dups, "rt")
    for d in sameseqs:
        logging.warn(d)
        g = re.match(r"\d+\t(.*)", d)
        logging.warn("%s", g.group(1))
        sp = g.group(1).split(",")
        for n in sp[1:]:
            print("<%s> <http://biohackathon.org/bh20-seq-schema/has_duplicate_sequence> <%s> ." % (n.strip(), sp[0].strip()))

orig = open(originalLabels, "rt")
print(orig.read())
