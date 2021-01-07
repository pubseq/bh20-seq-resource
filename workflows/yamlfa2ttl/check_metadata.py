import sys

import schema_salad.schema
import schema_salad.jsonld_context

from pyshex.evaluate import evaluate

path_yaml = sys.argv[1]
path_schema_yaml = sys.argv[2]
path_shex_rdf = sys.argv[3]

with open(path_schema_yaml, "rb") as f:
    cache = {
        "https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-schema.yml": f.read().decode("utf-8")
    }

metadata_schema = schema_salad.schema.load_schema(
    "https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-schema.yml",
    cache=cache
)

(document_loader, avsc_names, schema_metadata, metaschema_loader) = metadata_schema

if not isinstance(avsc_names, schema_salad.avro.schema.Names):
    raise Exception(avsc_names)

with open(path_shex_rdf, "rb") as f:
    shex = f.read().decode("utf-8")

doc, metadata = schema_salad.schema.load_and_validate(document_loader, avsc_names, path_yaml, True)
g = schema_salad.jsonld_context.makerdf("workflow", doc, document_loader.ctx)
rslt, reason = evaluate(g, shex, doc["id"], "https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-shex.rdf#submissionShape")

# As part of QC make sure serialization works too, this will raise
# an exception if there are invalid URIs.
g.serialize(format="ntriples")

if not rslt:
    raise Exception(reason)
