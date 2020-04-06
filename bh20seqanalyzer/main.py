import argparse
import arvados
import time
import subprocess
import tempfile
import json

def start_analysis(api, collection, analysis_project, workflow_uuid):
    project = api.groups().create(body={
        "group_class": "project",
        "name": "Analysis of %s" % collection["name"],
        "owner_uuid": analysis_project,
    }, ensure_unique_name=True).execute()

    with tempfile.NamedTemporaryFile() as tmp:
        inputobj = json.dumps({
            "sequence": {
                "class": "File",
                "location": "keep:%s/sequence.fasta" % collection["portable_data_hash"]
            },
            "metadata": {
                "class": "File",
                "location": "keep:%s/metadata.jsonld" % collection["portable_data_hash"]
            }
        }, indent=2)
        tmp.write(inputobj.encode('utf-8'))
        tmp.flush()
        cmd = ["arvados-cwl-runner",
               "--submit",
               "--no-wait",
               "--debug",
               "--project-uuid=%s" % project["uuid"],
               "arvwf:%s" % workflow_uuid,
               tmp.name]
        print("Running %s" % ' '.join(cmd))
        comp = subprocess.run(cmd, capture_output=True)
    if comp.returncode != 0:
        print(comp.stderr.decode('utf-8'))
    else:
        api.collections().update(uuid=collection["uuid"], body={"owner_uuid": project['uuid']}).execute()

def main():
    parser = argparse.ArgumentParser(description='Analyze collections uploaded to a project')
    parser.add_argument('--uploader-project', type=str, default='lugli-j7d0g-n5clictpuvwk8aa', help='')
    parser.add_argument('--analysis-project', type=str, default='lugli-j7d0g-y4k4uswcqi3ku56', help='')
    parser.add_argument('--workflow-uuid', type=str, default='lugli-7fd4e-mqfu9y3ofnpnho1', help='')
    args = parser.parse_args()

    api = arvados.api()

    while True:
        new_collections = api.collections().list(filters=[['owner_uuid', '=', args.uploader_project]]).execute()
        for c in new_collections["items"]:
            start_analysis(api, c, args.analysis_project, args.workflow_uuid)
        time.sleep(10)
