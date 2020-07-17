import arvados
import arvados.util
import arvados.keep
import ruamel.yaml

api = arvados.api()
keepclient = arvados.keep.KeepClient(api_client=api)

UPLOADER_PROJECT = 'lugli-j7d0g-n5clictpuvwk8aa'
VALIDATED_PROJECT = 'lugli-j7d0g-5ct8p1i1wrgyjvp'

delete_patterns = [
    "%missing%`collection_location`%",
    "%missing%`technology`%",
    "%missing%`host_species`%",
    "%QC fail: alignment%",
    "%does not look like a valid URI%",
    "%Duplicate of%",
    "%No matching triples found for predicate obo:NCIT_C42781%",
    "%does not look like a valid URI%"
    ]

revalidate_patterns = [
    "%missing%`license`%",
    "%QC fail%"
]

for p in delete_patterns:
    c = arvados.util.list_all(api.collections().list, filters=[
        ["owner_uuid", "=", UPLOADER_PROJECT],
        ["properties.errors", "like", p]])
    for i in c:
        print("trashing %s %s" % (i["uuid"], i["properties"].get("sequence_label")))
        api.collections().delete(uuid=i["uuid"]).execute()

for p in revalidate_patterns:
    c = arvados.util.list_all(api.collections().list, filters=[
        ["owner_uuid", "=", UPLOADER_PROJECT],
        ["properties.errors", "like", p]])
    for i in c:
        print("clearing status %s %s" % (i["uuid"], i["properties"].get("sequence_label")))
        pr = i["properties"]
        if "status" in pr:
            del pr["status"]
        if "errors" in pr:
            del pr["errors"]
        api.collections().update(uuid=i["uuid"], body={"properties": pr}).execute()

c = arvados.util.list_all(api.collections().list, filters=[
    ["owner_uuid", "=", VALIDATED_PROJECT],
    ["properties.sequence_label", "exists", False]])
for i in c:
    col = arvados.collection.Collection(i["uuid"], api_client=api, keep_client=keepclient)
    with col.open("metadata.yaml") as md:
        metadata_content = ruamel.yaml.round_trip_load(md)
    colprop = col.get_properties()
    colprop["sequence_label"] = metadata_content["sample"]["sample_id"]

    print("fixing sequence label %s %s" % (i["uuid"], colprop.get("sequence_label")))
    api.collections().update(uuid=i["uuid"], body={"properties": colprop}).execute()
