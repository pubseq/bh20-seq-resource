import sys
import arvados
import arvados.collection

from datetime import datetime

date_time_str = '2020-08-20'
date_time_obj = datetime.strptime(date_time_str, '%Y-%m-%d')

api = arvados.api()
keepclient = arvados.keep.KeepClient(api_client=api)

validated = arvados.util.list_all(api.collections().list, filters=[
    ["owner_uuid", "=", sys.argv[1]],
#    ["properties.status", "=", "validated"]
])

# validated.sort(key=lambda v: v["portable_data_hash"])

num_sample_deleted = 0
for item in validated:
    sequence_label = item['properties']["sequence_label"]

    # The SRA samples start with SRR or ERR
    if not sequence_label.startswith('SRR') and not sequence_label.startswith('ERR'):
        created_at_obj = datetime.strptime(item["created_at"], '%Y-%m-%dT%H:%M:%S.%fZ')
        # print(item, created_at_obj)

        if created_at_obj < date_time_obj:
            api.collections().delete(uuid=item['current_version_uuid']).execute()
            num_sample_deleted += 1
            print(sequence_label)

print('num_sample_deleted: {}'.format(num_sample_deleted))
