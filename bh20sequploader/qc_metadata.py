import schema_salad.schema
import schema_salad.ref_resolver
import schema_salad.jsonld_context
import logging
import pkg_resources
import logging
import traceback
from rdflib import Graph, Namespace
from pyshex.evaluate import evaluate


def qc_metadata(metadatafile):
    log = logging.getLogger(__name__ )
    schema_resource = pkg_resources.resource_stream(__name__, "bh20seq-schema.yml")
    cache = {"https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-schema.yml": schema_resource.read().decode("utf-8")}
    (document_loader,
     avsc_names,
     schema_metadata,
     metaschema_loader) = schema_salad.schema.load_schema("https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-schema.yml", cache=cache)

    shex = pkg_resources.resource_stream(__name__, "bh20seq-shex.rdf").read().decode("utf-8")

    if not isinstance(avsc_names, schema_salad.avro.schema.Names):
        print(avsc_names)
        return False

    try:
        doc, metadata = schema_salad.schema.load_and_validate(document_loader, avsc_names, metadatafile, True)
        g = schema_salad.jsonld_context.makerdf("workflow", doc, document_loader.ctx)
        rslt, reason = evaluate(g, shex, doc["id"], "https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-shex.rdf#submissionShape")

        if not rslt:
            log.debug(reason)
            print(reason)

        return rslt
    except Exception as e:
        traceback.print_exc()
        log.warn(e)
    return False
