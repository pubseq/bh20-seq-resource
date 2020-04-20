import schema_salad.schema
import schema_salad.jsonld_context

metadataSchema = '$(inputs.metadataSchema.path)'
metadata = $(inputs.metadata)
subjects = $(inputs.subjects)

(document_loader,
 avsc_names,
 schema_metadata,
 metaschema_loader) = schema_salad.schema.load_schema(metadataSchema)

for i, m in enumerate(metadata):
    doc, metadata = schema_salad.schema.load_and_validate(document_loader, avsc_names, m["path"], True)
    doc["id"] = subjects[i]
    g = schema_salad.jsonld_context.makerdf(subjects[i], doc, document_loader.ctx)
    print(g.serialize(format="ntriples").decode("utf-8"))
