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
from bh20sequploader.qc_fasta import qc_fasta
import pkg_resources
from schema_salad.sourceline import add_lc_filename

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
        try:
            metadata_content = ruamel.yaml.round_trip_load(col.open("metadata.yaml"))
            metadata_content["id"] = "http://arvados.org/keep:%s/metadata.yaml" % collection["portable_data_hash"]
            sample_id = metadata_content["sample"]["sample_id"]
            add_lc_filename(metadata_content, metadata_content["id"])
            valid = qc_metadata(metadata_content) and valid
        except Exception as e:
            logging.warn(e)
            valid = False
        if not valid:
            logging.warn("Failed metadata qc")

    if valid:
        try:
            tgt = None
            paired = {"reads_1.fastq": "reads.fastq", "reads_1.fastq.gz": "reads.fastq.gz"}
            for n in ("sequence.fasta", "reads.fastq", "reads.fastq.gz", "reads_1.fastq", "reads_1.fastq.gz"):
                if n not in col:
                    continue
                with col.open(n, 'rb') as qf:
                    tgt = qc_fasta(qf)[0]
                    if tgt != n and tgt != paired.get(n):
                        logging.info("Expected %s but magic says it should be %s", n, tgt)
                        valid = False
                    elif tgt in ("reads.fastq", "reads.fastq.gz", "reads_1.fastq", "reads_1.fastq.gz"):
                        start_fastq_to_fasta(api, collection, fastq_project, fastq_workflow_uuid, n, sample_id)
                        return False
            if tgt is None:
                valid = False
                logging.warn("Upload '%s' does not contain sequence.fasta, reads.fastq or reads_1.fastq", collection["name"])
        except ValueError as v:
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
        api.collections().update(uuid=collection["uuid"], body={
            "owner_uuid": validated_project,
            "name": "%s (%s)" % (collection["name"], time.asctime(time.gmtime()))}).execute()
    else:
        # It is invalid, delete it.
        logging.warn("Suggest deleting '%s'" % collection["name"])
        #api.collections().delete(uuid=collection["uuid"]).execute()

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
               "--project-uuid=%s" % project["uuid"],
               "arvwf:%s" % workflow_uuid,
               tmp.name]
        logging.info("Running %s" % ' '.join(cmd))
        comp = subprocess.run(cmd, capture_output=True)
    logging.info("Submitted %s", comp.stdout)
    if comp.returncode != 0:
        logging.error(comp.stderr.decode('utf-8'))

    return project


def start_fastq_to_fasta(api, collection,
                         analysis_project,
                         fastq_workflow_uuid,
                         tgt,
                         sample_id):

    params = {
        "metadata": {
            "class": "File",
            "location": "keep:%s/metadata.yaml" % collection["portable_data_hash"]
        },
        "ref_fasta": {
            "class": "File",
            "location": "keep:ffef6a3b77e5e04f8f62a7b6f67264d1+556/SARS-CoV2-NC_045512.2.fasta"
        },
        "sample_id": sample_id
    }

    if tgt.startswith("reads.fastq"):
        params["fastq_forward"] = {
            "class": "File",
            "location": "keep:%s/%s" % (collection["portable_data_hash"], tgt)
        }
    elif tgt.startswith("reads_1.fastq"):
        params["fastq_forward"] = {
            "class": "File",
            "location": "keep:%s/reads_1.%s" % (collection["portable_data_hash"], tgt[8:])
        }
        params["fastq_reverse"] = {
            "class": "File",
            "location": "keep:%s/reads_2.%s" % (collection["portable_data_hash"], tgt[8:])
        }

    newproject = run_workflow(api, analysis_project, fastq_workflow_uuid, "FASTQ to FASTA", params)
    api.collections().update(uuid=collection["uuid"],
                             body={"owner_uuid": newproject["uuid"]}).execute()

def start_pangenome_analysis(api,
                             analysis_project,
                             pangenome_workflow_uuid,
                             validated_project,
                             schema_ref,
                             exclude_list):
    validated = arvados.util.list_all(api.collections().list, filters=[["owner_uuid", "=", validated_project]])
    inputobj = {
        "inputReads": [],
        "metadata": [],
        "subjects": [],
        "metadataSchema": {
            "class": "File",
            "location": schema_ref
        },
        "exclude": {
            "class": "File",
            "location": exclude_list
        }
    }
    validated.sort(key=lambda v: v["portable_data_hash"])
    for v in validated:
        inputobj["inputReads"].append({
            "class": "File",
            "location": "keep:%s/sequence.fasta" % v["portable_data_hash"]
        })
        inputobj["metadata"].append({
            "class": "File",
            "location": "keep:%s/metadata.yaml" % v["portable_data_hash"]
        })
        inputobj["subjects"].append("http://collections.lugli.arvadosapi.com/c=%s/sequence.fasta" % v["portable_data_hash"])
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


def upload_schema(api, workflow_def_project):
    schema_resource = pkg_resources.resource_stream('bh20sequploader.qc_metadata', "bh20seq-schema.yml")
    c = arvados.collection.Collection()
    with c.open("schema.yml", "wb") as f:
        f.write(schema_resource.read())
    pdh = c.portable_data_hash()
    wd = api.collections().list(filters=[["owner_uuid", "=", workflow_def_project],
                                         ["portable_data_hash", "=", pdh]]).execute()
    if len(wd["items"]) == 0:
        c.save_new(owner_uuid=workflow_def_project, name="Metadata schema", ensure_unique_name=True)
    return "keep:%s/schema.yml" % pdh


def main():
    parser = argparse.ArgumentParser(description='Analyze collections uploaded to a project')
    parser.add_argument('--uploader-project', type=str, default='lugli-j7d0g-n5clictpuvwk8aa', help='')
    parser.add_argument('--pangenome-analysis-project', type=str, default='lugli-j7d0g-y4k4uswcqi3ku56', help='')
    parser.add_argument('--fastq-project', type=str, default='lugli-j7d0g-xcjxp4oox2u1w8u', help='')
    parser.add_argument('--validated-project', type=str, default='lugli-j7d0g-5ct8p1i1wrgyjvp', help='')
    parser.add_argument('--workflow-def-project', type=str, default='lugli-j7d0g-5hswinmpyho8dju', help='')

    parser.add_argument('--pangenome-workflow-uuid', type=str, default='lugli-7fd4e-mqfu9y3ofnpnho1', help='')
    parser.add_argument('--fastq-workflow-uuid', type=str, default='lugli-7fd4e-2zp9q4jo5xpif9y', help='')

    parser.add_argument('--exclude-list', type=str, default='keep:lugli-4zz18-tzzhcm6hrf8ci8d/exclude.txt', help='')

    parser.add_argument('--latest-result-collection', type=str, default='lugli-4zz18-z513nlpqm03hpca', help='')
    parser.add_argument('--kickoff', action="store_true")
    parser.add_argument('--once', action="store_true")
    args = parser.parse_args()

    api = arvados.api()



    schema_ref = upload_schema(api, args.workflow_def_project)

    if args.kickoff:
        logging.info("Starting a single analysis run")
        start_pangenome_analysis(api,
                                 args.pangenome_analysis_project,
                                 args.pangenome_workflow_uuid,
                                 args.validated_project,
                                 schema_ref,
                                 args.exclude_list)
        return

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
                                     args.validated_project,
                                     schema_ref,
                                     args.exclude_list)

        copy_most_recent_result(api,
                                args.pangenome_analysis_project,
                                args.latest_result_collection)

        if args.once:
            break
        time.sleep(15)
