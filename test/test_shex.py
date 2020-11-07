import schema_salad.schema
import schema_salad.ref_resolver
import schema_salad.jsonld_context
import pkg_resources
import logging
import traceback
# from rdflib import Graph, Namespace
from pyshex.evaluate import evaluate
import unittest

class TestStringMethods(unittest.TestCase):

    def test_schema(self):
        with open("bh20sequploader/bh20seq-schema.yml") as schema_resource:
            metadata_schema = schema_salad.schema.load_schema("bh20sequploader/bh20seq-schema.yml")
            (document_loader,
             avsc_names,
             schema_metadata,
             metaschema_loader) = metadata_schema
            print(metadata_schema)
            assert(isinstance(avsc_names, schema_salad.avro.schema.Names))
            metadatafile = "test/data/input/TN_UT2.yaml"
            doc, metadata = schema_salad.schema.load_and_validate(document_loader, avsc_names, metadatafile, True)
            print(doc)
            g = schema_salad.jsonld_context.makerdf("workflow", doc, document_loader.ctx)
            shex = pkg_resources.resource_stream(__name__, "../bh20sequploader/bh20seq-shex.rdf").read().decode("utf-8")
            # Note the https link simply acts as a URI descriptor (it does not fetch)
            rslt, reason = evaluate(g, shex, doc["id"], "https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-shex.rdf#submissionShape")

            g.serialize(format="ntriples")

            if not rslt:
                raise Exception(reason)

if __name__ == '__main__':
    unittest.main()
