import arvados
import arvados.util

api = arvados.api()

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
        ["owner_uuid", "=", "lugli-j7d0g-n5clictpuvwk8aa"],
        ["properties.errors", "like", p]])
    for i in c:
        print("trashing %s %s" % (i["uuid"], i["properties"].get("sequence_label")))
        api.collections().delete(uuid=i["uuid"]).execute()

for p in revalidate_patterns:
    c = arvados.util.list_all(api.collections().list, filters=[
        ["owner_uuid", "=", "lugli-j7d0g-n5clictpuvwk8aa"],
        ["properties.errors", "like", p]])
    for i in c:
        print("clearing status %s %s" % (i["uuid"], i["properties"].get("sequence_label")))
        pr = i["properties"]
        if "status" in pr:
            del pr["status"]
        if "errors" in pr:
            del pr["errors"]
        api.collections().update(uuid=i["uuid"], body={"properties": pr}).execute()
