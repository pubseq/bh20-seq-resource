import argparse
import arvados
import arvados.collection
import time
import subprocess
import tempfile
import json
import logging
import ruamel.yaml
from bh20sequploader.qc_metadata import qc_metadata

logging.basicConfig(format="[%(asctime)s] %(levelname)s %(message)s", datefmt="%Y-%m-%d %H:%M:%S",
                    level=logging.INFO)
logging.getLogger("googleapiclient.discovery").setLevel(logging.WARN)

def validate_upload(api, collection, validated_project,
                    fastq_project, fastq_workflow_uuid):
    col = arvados.collection.Collection(collection["uuid"])

    # validate the collection here.  Check metadata, etc.
    valid = True

    if "metadata.yaml" not in col:
        logging.warn("Upload '%s' missing metadata.yaml", collection["name"])
        valid = False
    else:
        metadata_content = ruamel.yaml.round_trip_load(col.open("metadata.yaml"))
        #valid = qc_metadata(metadata_content) and valid
        if not valid:
            logging.warn("Failed metadata qc")

    if valid:
        if "sequence.fasta" not in col:
            if "reads.fastq" in col:
                start_fastq_to_fasta(api, collection, fastq_project, fastq_workflow_uuid)
                return False
            else:
                valid = False
                logging.warn("Upload '%s' missing sequence.fasta", collection["name"])

    dup = api.collections().list(filters=[["owner_uuid", "=", validated_project],
                                          ["portable_data_hash", "=", col.portable_data_hash()]]).execute()
    if dup["items"]:
        # This exact collection has been uploaded before.
        valid = False
        logging.warn("Upload '%s' is duplicate" % collection["name"])

    if valid:
        logging.info("Added '%s' to validated sequences" % collection["name"])
        # Move it to the "validated" project to be included in the next analysis
        api.collections().update(uuid=collection["uuid"], body={
            "owner_uuid": validated_project,
            "name": "%s (%s)" % (collection["name"], time.asctime(time.gmtime()))}).execute()
    else:
        # It is invalid, delete it.
        logging.warn("Deleting '%s'" % collection["name"])
        api.collections().delete(uuid=collection["uuid"]).execute()

    return valid


def run_workflow(api, parent_project, workflow_uuid, name, inputobj):
    project = api.groups().create(body={
        "group_class": "project",
        "name": name,
        "owner_uuid": parent_project,
    }, ensure_unique_name=True).execute()

    with tempfile.NamedTemporaryFile() as tmp:
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

    return project


def start_fastq_to_fasta(api, collection,
                         analysis_project,
                         fastq_workflow_uuid):
    newproject = run_workflow(api, analysis_project, fastq_workflow_uuid, "FASTQ to FASTA", {
        "fastq_forward": {
            "class": "File",
            "location": "keep:%s/reads.fastq" % collection["portable_data_hash"]
        },
        "metadata": {
            "class": "File",
            "location": "keep:%s/metadata.yaml" % collection["portable_data_hash"]
        },
        "ref_fasta": {
            "class": "File",
            "location": "keep:ffef6a3b77e5e04f8f62a7b6f67264d1+556/SARS-CoV2-NC_045512.2.fasta"
        }
    })
    api.collections().update(uuid=collection["uuid"],
                             body={"owner_uuid": newproject["uuid"]}).execute()

def start_pangenome_analysis(api,
                             analysis_project,
                             pangenome_workflow_uuid,
                             validated_project):
    validated = arvados.util.list_all(api.collections().list, filters=[["owner_uuid", "=", validated_project]])
    inputobj = {
        "inputReads": []
    }
    for v in validated:
        inputobj["inputReads"].append({
            "class": "File",
            "location": "keep:%s/sequence.fasta" % v["portable_data_hash"]
        })
    run_workflow(api, analysis_project, pangenome_workflow_uuid, "Pangenome analysis", inputobj)


def get_workflow_output_from_project(api, uuid):
    cr = api.container_requests().list(filters=[['owner_uuid', '=', uuid],
                                                ["requesting_container_uuid", "=", None]]).execute()
    if cr["items"] and cr["items"][0]["output_uuid"]:
        return cr["items"][0]
    else:
        return None


def copy_most_recent_result(api, analysis_project, latest_result_uuid):
    most_recent_analysis = api.groups().list(filters=[['owner_uuid', '=', analysis_project]],
                                                  order="created_at desc", limit=1).execute()
    for m in most_recent_analysis["items"]:
        wf = get_workflow_output_from_project(api, m["uuid"])
        if wf:
            src = api.collections().get(uuid=wf["output_uuid"]).execute()
            dst = api.collections().get(uuid=latest_result_uuid).execute()
            if src["portable_data_hash"] != dst["portable_data_hash"]:
                logging.info("Copying latest result from '%s' to %s", m["name"], latest_result_uuid)
                api.collections().update(uuid=latest_result_uuid,
                                         body={"manifest_text": src["manifest_text"],
                                               "description": "Result from %s %s" % (m["name"], wf["uuid"])}).execute()
            break


def move_fastq_to_fasta_results(api, analysis_project, uploader_project):
    projects = api.groups().list(filters=[['owner_uuid', '=', analysis_project],
                                          ["properties.moved_output", "!=", True]],
                                 order="created_at desc",).execute()
    for p in projects["items"]:
        wf = get_workflow_output_from_project(api, p["uuid"])
        if wf:
            logging.info("Moving completed fastq2fasta result %s back to uploader project", wf["output_uuid"])
            api.collections().update(uuid=wf["output_uuid"],
                                     body={"owner_uuid": uploader_project}).execute()
            p["properties"]["moved_output"] = True
            api.groups().update(uuid=p["uuid"], body={"properties": p["properties"]}).execute()


def main():
    parser = argparse.ArgumentParser(description='Analyze collections uploaded to a project')
    parser.add_argument('--uploader-project', type=str, default='lugli-j7d0g-n5clictpuvwk8aa', help='')
    parser.add_argument('--pangenome-analysis-project', type=str, default='lugli-j7d0g-y4k4uswcqi3ku56', help='')
    parser.add_argument('--fastq-project', type=str, default='lugli-j7d0g-xcjxp4oox2u1w8u', help='')
    parser.add_argument('--validated-project', type=str, default='lugli-j7d0g-5ct8p1i1wrgyjvp', help='')

    parser.add_argument('--pangenome-workflow-uuid', type=str, default='lugli-7fd4e-mqfu9y3ofnpnho1', help='')
    parser.add_argument('--fastq-workflow-uuid', type=str, default='lugli-7fd4e-2zp9q4jo5xpif9y', help='')

    parser.add_argument('--latest-result-collection', type=str, default='lugli-4zz18-z513nlpqm03hpca', help='')
    args = parser.parse_args()

    api = arvados.api()

    logging.info("Starting up, monitoring %s for uploads" % (args.uploader_project))

    while True:
        move_fastq_to_fasta_results(api, args.fastq_project, args.uploader_project)

        new_collections = api.collections().list(filters=[['owner_uuid', '=', args.uploader_project]]).execute()
        at_least_one_new_valid_seq = False
        for c in new_collections["items"]:
            at_least_one_new_valid_seq = validate_upload(api, c,
                                                         args.validated_project,
                                                         args.fastq_project,
                                                         args.fastq_workflow_uuid) or at_least_one_new_valid_seq

        if at_least_one_new_valid_seq:
            start_pangenome_analysis(api,
                                     args.pangenome_analysis_project,
                                     args.pangenome_workflow_uuid,
                                     args.validated_project)

        copy_most_recent_result(api,
                                args.pangenome_analysis_project,
                                args.latest_result_collection)

        time.sleep(15)
