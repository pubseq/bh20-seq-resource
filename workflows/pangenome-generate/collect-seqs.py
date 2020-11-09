import sys
import arvados
import json
import shutil
import arvados.collection
import ruamel.yaml
import schema_salad.schema
import schema_salad.jsonld_context
from schema_salad.sourceline import add_lc_filename

api = arvados.api()
keepclient = arvados.keep.KeepClient(api_client=api)

validated = arvados.util.list_all(api.collections().list, filters=[
    ["owner_uuid", "=", sys.argv[1]],
    ["properties.status", "=", "validated"]])

validated.sort(key=lambda v: v["portable_data_hash"])

relabeled_fasta = open("relabeledSeqs.fasta", "wt")
merged_metadata = open("mergedMetadata.ttl", "wt")

metadataSchema = sys.argv[2]

blacklist = set()
if len(sys.argv) > 3:
    with open(sys.argv[3]) as bl:
        for l in bl:
            blacklist.add(l.strip())

(document_loader,
 avsc_names,
 schema_metadata,
 metaschema_loader) = schema_salad.schema.load_schema(metadataSchema)


for item in validated:
    pdh = item["portable_data_hash"]
    with arvados.collection.CollectionReader(pdh, api_client=api, keep_client=keepclient) as col:
        with col.open("sequence.fasta", "rt") as fa:
            subject = "http://covid19.genenetwork.org/resource/%s" % pdh
            label = fa.readline().strip()
            merged_metadata.write("<%s> <http://biohackathon.org/bh20-seq-schema/original_fasta_label> \"%s\" .\n" % (subject, label[1:].replace('"', '\\"')))
            skip = (subject in blacklist or label[1:] in blacklist)
            if skip:
                merged_metadata.write("<%s> <http://biohackathon.org/bh20-seq-schema/excluded_from_graph> \"true\"^^<http://www.w3.org/2001/XMLSchema#boolean> .\n" % subject)
            if not skip:
                relabeled_fasta.write(">"+subject+"\n")
            data = fa.read(8096)
            while data:
                if not skip:
                    relabeled_fasta.write(data)
                endswithnewline = data.endswith("\n")
                data = fa.read(8096)
            if not skip and not endswithnewline:
                relabeled_fasta.write("\n")

        with col.open("metadata.yaml", "rt") as md:
            metadata_content = ruamel.yaml.round_trip_load(md)
        metadata_content["id"] = subject
        add_lc_filename(metadata_content, metadata_content["id"])
        doc, metadata = schema_salad.schema.load_and_validate(document_loader, avsc_names, metadata_content, False, False)
        g = schema_salad.jsonld_context.makerdf(subject, doc, document_loader.ctx)
        merged_metadata.write(g.serialize(format="ntriples").decode("utf-8"))


shutil.rmtree(".cache")
