import yamale

## NOTE: this is just a DUMMY. Everything about this can and will change
def qc_metadata(metadatafile):
    print("Start metadata validation...")
    schema = yamale.make_schema('../example/dummyschema.yaml')
    data = yamale.make_data(metadatafile)
    # Validate data against the schema. Throws a ValueError if data is invalid.
    yamale.validate(schema, data)
    print("...complete!")

#qc_metadata("../example/metadata.yaml")

