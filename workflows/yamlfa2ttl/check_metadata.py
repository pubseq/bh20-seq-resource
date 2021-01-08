import sys

import schema_salad.schema
import schema_salad.jsonld_context

from pyshex.evaluate import evaluate

path_yaml = sys.argv[1]
path_fasta = sys.argv[2]
path_schema_yaml = sys.argv[3]
path_shex_rdf = sys.argv[4]

with open(path_schema_yaml, "rb") as f:
    cache = {
        "https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-schema.yml": f.read().decode(
            "utf-8")
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
rslt, reason = evaluate(
    g, shex, doc["id"],
    "https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-shex.rdf#submissionShape"
)

# As part of QC make sure serialization works too, this will raise
# an exception if there are invalid URIs.
g.serialize(format="ntriples")

if not rslt:
    raise Exception(reason)

# The sample_id in the FASTA header has to equal to the sample_id in the YAML file
sample_id_from_metadata = metadata['sample']['sample_id']

sample_id_from_fasta = ''

with open(path_fasta) as f:
    for line in f:
        sample_id_from_fasta = line.strip().split(' ')[0][1:]
        break

if sample_id_from_metadata != sample_id_from_fasta:
    raise ValueError(
        f"sample_id in the YAML file '{sample_id_from_metadata}' is different from the sample_id in the FASTA '{sample_id_from_fasta}'"
    )
