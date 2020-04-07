import argparse
import arvados
import arvados.collection
import time
import subprocess
import tempfile
import json
import logging

logging.basicConfig(format="[%(asctime)s] %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S",
                    level=logging.INFO)
logging.getLogger("googleapiclient.discovery").setLevel(logging.WARN)

def validate_upload(api, collection, validated_project, latest_result_uuid):
    col = arvados.collection.Collection(collection["uuid"])

    # validate the collection here.  Check metadata, etc.
    valid = True

    if "sequence.fasta" not in col:
        valid = False
        logging.warn("Upload '%s' missing sequence.fasta", collection["name"])
    if "metadata.jsonld" not in col:
        logging.warn("Upload '%s' missing metadata.jsonld", collection["name"])
        valid = False

    dup = api.collections().list(filters=[["owner_uuid", "=", validated_project],
                                          ["portable_data_hash", "=", col.portable_data_hash()]]).execute()
    if dup["items"]:
        # This exact collection has been uploaded before.
        valid = False
        logging.warn("Upload '%s' is duplicate" % collection["name"])

    if valid:
        logging.info("Added '%s' to validated sequences" % collection["name"])
        # Move it to the "validated" project to be included in the next analysis
        api.collections().update(uuid=collection["uuid"], body={"owner_uuid": validated_project}).execute()
    else:
        # It is invalid, delete it.
        logging.warn("Deleting '%s'" % collection["name"])
        api.collections().delete(uuid=collection["uuid"]).execute()

    return valid

def start_analysis(api,
                   analysis_project,
                   workflow_uuid,
                   validated_project):

    project = api.groups().create(body={
        "group_class": "project",
        "name": "Pangenome analysis",
        "owner_uuid": analysis_project,
    }, ensure_unique_name=True).execute()

    validated = arvados.util.list_all(api.collections().list, filters=[["owner_uuid", "=", validated_project]])

    with tempfile.NamedTemporaryFile() as tmp:
        inputobj = {
            "inputReads": []
        }
        for v in validated:
            inputobj["inputReads"].append({
                "class": "File",
                "location": "keep:%s/sequence.fasta" % v["portable_data_hash"]
            })
        tmp.write(json.dumps(inputobj, indent=2).encode('utf-8'))
        tmp.flush()
        cmd = ["arvados-cwl-runner",
               "--submit",
               "--no-wait",
               "--debug",
               "--project-uuid=%s" % project["uuid"],
               "arvwf:%s" % workflow_uuid,
               tmp.name]
        logging.info("Running %s" % ' '.join(cmd))
        comp = subprocess.run(cmd, capture_output=True)
    if comp.returncode != 0:
        logging.error(comp.stderr.decode('utf-8'))


def copy_most_recent_result(api, analysis_project, latest_result_uuid):
    most_recent_analysis = api.groups().list(filters=[['owner_uuid', '=', analysis_project]],
                                                  order="created_at desc", limit=1).execute()
    for m in most_recent_analysis["items"]:
        cr = api.container_requests().list(filters=[['owner_uuid', '=', m["uuid"]],
                                                    ["requesting_container_uuid", "=", None]]).execute()
        if cr["items"] and cr["items"][0]["output_uuid"]:
            wf = cr["items"][0]
            src = api.collections().get(uuid=wf["output_uuid"]).execute()
            dst = api.collections().get(uuid=latest_result_uuid).execute()
            if src["portable_data_hash"] != dst["portable_data_hash"]:
                logging.info("Copying latest result from '%s' to %s", m["name"], latest_result_uuid)
                api.collections().update(uuid=latest_result_uuid,
                                         body={"manifest_text": src["manifest_text"],
                                               "description": "latest result from %s %s" % (m["name"], wf["uuid"])}).execute()
            break


def main():
    parser = argparse.ArgumentParser(description='Analyze collections uploaded to a project')
    parser.add_argument('--uploader-project', type=str, default='lugli-j7d0g-n5clictpuvwk8aa', help='')
    parser.add_argument('--analysis-project', type=str, default='lugli-j7d0g-y4k4uswcqi3ku56', help='')
    parser.add_argument('--validated-project', type=str, default='lugli-j7d0g-5ct8p1i1wrgyjvp', help='')
    parser.add_argument('--workflow-uuid', type=str, default='lugli-7fd4e-mqfu9y3ofnpnho1', help='')
    parser.add_argument('--latest-result-uuid', type=str, default='lugli-4zz18-z513nlpqm03hpca', help='')
    args = parser.parse_args()

    api = arvados.api()

    logging.info("Starting up, monitoring %s for uploads" % (args.uploader_project))

    while True:
        new_collections = api.collections().list(filters=[['owner_uuid', '=', args.uploader_project]]).execute()
        at_least_one_new_valid_seq = False
        for c in new_collections["items"]:
            at_least_one_new_valid_seq = validate_upload(api, c, args.validated_project) or at_least_one_new_valid_seq

        if at_least_one_new_valid_seq:
            start_analysis(api, args.analysis_project,
                           args.workflow_uuid,
                           args.validated_project)

        copy_most_recent_result(api, args.analysis_project, args.latest_result_uuid)

        time.sleep(10)
