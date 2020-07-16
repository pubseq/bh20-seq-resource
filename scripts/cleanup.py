import arvados
import arvados.util

api = arvados.api()

patterns = [
    "%missing%`collection_location`%",
    "%missing%`technology`%",
    "%missing%`host_species`%",
    "%QC fail: alignment%",
    "%does not look like a valid URI%",
    ]

for p in patterns:
    c = arvados.util.list_all(api.collections().list, filters=[
        ["owner_uuid", "=", "lugli-j7d0g-n5clictpuvwk8aa"],
        ["properties.errors", "like", p]])
    for i in c:
        print("trashing %s %s" % (i["uuid"], i["properties"].get("sequence_label")))
        api.collections().delete(uuid=i["uuid"]).execute()
