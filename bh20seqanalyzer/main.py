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

class SeqAnalyzer:

    def __init__(self, api, keepclient,
                 uploader_project,
                 pangenome_analysis_project,
                 fastq_project,
                 validated_project,
                 workflow_def_project,
                 pangenome_workflow_uuid,
                 fastq_workflow_uuid,
                 exclude_list,
                 latest_result_collection):
        self.api = api
        self.keepclient = keepclient
        self.uploader_project = uploader_project
        self.pangenome_analysis_project = pangenome_analysis_project
        self.fastq_project = fastq_project
        self.validated_project = validated_project
        self.workflow_def_project = workflow_def_project
        self.pangenome_workflow_uuid = pangenome_workflow_uuid
        self.fastq_workflow_uuid = fastq_workflow_uuid
        self.exclude_list = exclude_list
        self.latest_result_uuid = latest_result_collection
        self.schema_ref = None

    def validate_upload(self, collection, revalidate):
        col = arvados.collection.Collection(collection["uuid"], api_client=self.api, keep_client=self.keepclient)

        if not revalidate and collection["properties"].get("status") in ("validated", "rejected"):
            return False

        # validate the collection here.  Check metadata, etc.
        logging.info("Validating upload '%s' (%s)" % (collection["name"], collection["uuid"]))

        errors = []

        if collection["owner_uuid"] != self.validated_project:
            dup = self.api.collections().list(filters=[["owner_uuid", "=", self.validated_project],
                                                  ["portable_data_hash", "=", col.portable_data_hash()]]).execute()
            if dup["items"]:
                # This exact collection has been uploaded before.
                errors.append("Duplicate of %s" % ([d["uuid"] for d in dup["items"]]))

        if not errors:
            if "metadata.yaml" not in col:
                errors.append("Missing metadata.yaml", collection["name"])
            else:
                try:
                    with col.open("metadata.yaml") as md:
                        metadata_content = ruamel.yaml.round_trip_load(md)
                    metadata_content["id"] = "http://arvados.org/keep:%s/metadata.yaml" % collection["portable_data_hash"]
                    sample_id = metadata_content["sample"]["sample_id"]
                    add_lc_filename(metadata_content, metadata_content["id"])
                    valid = qc_metadata(metadata_content)
                    if not valid:
                        errors.append("Failed metadata qc")
                except Exception as e:
                    errors.append(str(e))

        if not errors:
            try:
                tgt = None
                paired = {"reads_1.fastq": "reads.fastq", "reads_1.fastq.gz": "reads.fastq.gz"}
                for n in ("sequence.fasta", "reads.fastq", "reads.fastq.gz", "reads_1.fastq", "reads_1.fastq.gz"):
                    if n not in col:
                        continue
                    with col.open(n, 'rb') as qf:
                        tgt, seqlabel, seq_type = qc_fasta(qf)
                        if tgt != n and tgt != paired.get(n):
                            errors.append("Expected %s but magic says it should be %s", n, tgt)
                        elif tgt in ("reads.fastq", "reads.fastq.gz", "reads_1.fastq", "reads_1.fastq.gz"):
                            self.start_fastq_to_fasta(collection, n, sample_id)
                            return False

                        # If it is a FASTA
                        if sample_id != seqlabel:
                            errors.append("Expected sample_id == seqlabel, but %s != %s" % (sample_id, seqlabel))
                if tgt is None:
                    errors.append("Upload '%s' does not contain sequence.fasta, reads.fastq or reads_1.fastq", collection["name"])
            except Exception as v:
                errors.append(str(v))


        if not errors:
            # Move it to the "validated" project to be included in the next analysis
            if "errors" in collection["properties"]:
                del collection["properties"]["errors"]
            collection["properties"]["status"] = "validated"
            self.api.collections().update(uuid=collection["uuid"], body={
                "owner_uuid": self.validated_project,
                "name": "%s (%s)" % (collection["name"], time.asctime(time.gmtime())),
                "properties": collection["properties"]}).execute()
            logging.info("Added '%s' to validated sequences" % collection["name"])
            return True
        else:
            # It is invalid
            logging.warn("'%s' (%s) has validation errors: %s" % (
                collection["name"], collection["uuid"], "\n".join(errors)))
            collection["properties"]["status"] = "rejected"
            collection["properties"]["errors"] = errors
            self.api.collections().update(uuid=collection["uuid"], body={"properties": collection["properties"]}).execute()
            return False


    def run_workflow(self, parent_project, workflow_uuid, name, inputobj):
        project = self.api.groups().create(body={
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


    def start_fastq_to_fasta(self, collection,
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

        newproject = self.run_workflow(self.fastq_project, self.fastq_workflow_uuid, "FASTQ to FASTA", params)
        self.api.collections().update(uuid=collection["uuid"],
                                 body={"owner_uuid": newproject["uuid"]}).execute()

    def start_pangenome_analysis(self):

        if self.schema_ref is None:
            self.upload_schema()

        inputobj = {
            "metadataSchema": {
                "class": "File",
                "location": self.schema_ref
            },
            "exclude": {
                "class": "File",
                "location": self.exclude_list
            },
            "src_project": self.validated_project
        }

        self.run_workflow(self.pangenome_analysis_project, self.pangenome_workflow_uuid, "Pangenome analysis", inputobj)


    def get_workflow_output_from_project(self, uuid, named):
        cr = self.api.container_requests().list(filters=[['owner_uuid', '=', uuid],
                                                         ["requesting_container_uuid", "=", None],
                                                         ["name", "=", named]]).execute()
        if cr["items"] and cr["items"][0]["output_uuid"]:
            container = self.api.containers().get(uuid=cr["items"][0]["container_uuid"]).execute()
            if container["state"] == "Complete" and container["exit_code"] == 0:
                return cr["items"][0]
        return None


    def copy_most_recent_result(self):
        most_recent_analysis = self.api.groups().list(filters=[['owner_uuid', '=', self.pangenome_analysis_project]],
                                                      order="created_at desc").execute()
        for m in most_recent_analysis["items"]:
            wf = self.get_workflow_output_from_project(m["uuid"], "collect-seqs.cwl")
            if wf is None:
                continue
            src = self.api.collections().get(uuid=wf["output_uuid"]).execute()
            dst = self.api.collections().get(uuid=self.latest_result_uuid).execute()
            if src["portable_data_hash"] != dst["portable_data_hash"]:
                logging.info("Copying latest result from '%s' to %s", m["name"], self.latest_result_uuid)
                self.api.collections().update(uuid=self.latest_result_uuid,
                                         body={"manifest_text": src["manifest_text"],
                                               "description": "Result from %s %s" % (m["name"], wf["uuid"])}).execute()
            break


    def move_fastq_to_fasta_results(self):
        projects = arvados.util.list_all(self.api.groups().list,
                                         filters=[['owner_uuid', '=', self.fastq_project],
                                                ["properties.moved_output", "!=", True]],
                                         order="created_at asc")
        for p in projects:
            wf = self.get_workflow_output_from_project(p["uuid"], "fastq2fasta.cwl")
            if not wf:
                continue

            logging.info("Moving completed fastq2fasta result %s back to uploader project", wf["output_uuid"])

            col = arvados.collection.Collection(wf["output_uuid"], api_client=self.api, keep_client=self.keepclient)
            with col.open("metadata.yaml") as md:
                metadata_content = ruamel.yaml.round_trip_load(md)

            colprop = col.get_properties()
            colprop["sequence_label"] = metadata_content["sample"]["sample_id"]
            self.api.collections().update(uuid=wf["output_uuid"],
                                     body={"owner_uuid": self.uploader_project,
                                           "properties": colprop}).execute()

            p["properties"]["moved_output"] = True
            self.api.groups().update(uuid=p["uuid"], body={"properties": p["properties"]}).execute()


    def upload_schema(self):
        schema_resource = pkg_resources.resource_stream('bh20sequploader.qc_metadata', "bh20seq-schema.yml")
        c = arvados.collection.Collection(api_client=self.api, keep_client=self.keepclient)
        with c.open("schema.yml", "wb") as f:
            f.write(schema_resource.read())
        pdh = c.portable_data_hash()
        wd = self.api.collections().list(filters=[["owner_uuid", "=", self.workflow_def_project],
                                             ["portable_data_hash", "=", pdh]]).execute()
        if len(wd["items"]) == 0:
            c.save_new(owner_uuid=self.workflow_def_project, name="Metadata schema", ensure_unique_name=True)
        self.schema_ref = "keep:%s/schema.yml" % pdh


    def print_status(self, fmt):
        pending = arvados.util.list_all(self.api.collections().list, filters=[["owner_uuid", "=", self.uploader_project]])
        out = []
        status = {}
        for p in pending:
            prop = p["properties"]
            out.append(prop)
            if "status" not in prop:
                prop["status"] = "pending"
            prop["created_at"] = p["created_at"]
            prop["uuid"] = p["uuid"]
            status[prop["status"]] = status.get(prop["status"], 0) + 1
        if fmt == "html":
            print(
    """
    <html>
    <body>
    """)
            print("<p>Total collections in upload project %s</p>" % len(out))
            print("<p>Status %s</p>" % status)
            print(
    """
    <table>
    <tr><th>Collection</th>
    <th>Sequence label</th>
    <th>Status</th>
    <th>Errors</th></tr>
    """)
            for r in out:
                print("<tr valign='top'>")
                print("<td><a href='https://workbench.lugli.arvadosapi.com/collections/%s'>%s</a></td>" % (r["uuid"], r["uuid"]))
                print("<td>%s</td>" % r["sequence_label"])
                print("<td>%s</td>" % r["status"])
                print("<td><pre>%s</pre></td>" % "\n".join(r.get("errors", [])))
                print("</tr>")
            print(
    """
    </table>
    </body>
    </html>
    """)
        else:
            print(json.dumps(out, indent=2))

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
    parser.add_argument('--no-start-analysis', action="store_true")
    parser.add_argument('--once', action="store_true")
    parser.add_argument('--print-status', type=str, default=None)
    parser.add_argument('--revalidate', action="store_true", default=None)
    args = parser.parse_args()

    api = arvados.api()
    keepclient = arvados.keep.KeepClient(api_client=api)

    seqanalyzer = SeqAnalyzer(api, keepclient,
                              args.uploader_project,
                              args.pangenome_analysis_project,
                              args.fastq_project,
                              args.validated_project,
                              args.workflow_def_project,
                              args.pangenome_workflow_uuid,
                              args.fastq_workflow_uuid,
                              args.exclude_list,
                              args.latest_result_collection)

    if args.kickoff:
        logging.info("Starting a single analysis run")
        seqanalyzer.start_pangenome_analysis()
        return

    if args.print_status:
        seqanalyzer.print_status(args.print_status)
        exit(0)

    logging.info("Starting up, monitoring %s for uploads" % (args.uploader_project))

    while True:
        try:
            seqanalyzer.move_fastq_to_fasta_results()

            new_collections = arvados.util.list_all(api.collections().list, filters=[["owner_uuid", "=", args.uploader_project]])
            at_least_one_new_valid_seq = False
            for c in new_collections:
                at_least_one_new_valid_seq = seqanalyzer.validate_upload(c, args.revalidate) or at_least_one_new_valid_seq

            if at_least_one_new_valid_seq and not args.no_start_analysis:
                seqanalyzer.start_pangenome_analysis()

            seqanalyzer.copy_most_recent_result()
        except Exception as e:
            logging.exeception("Error in main loop")

        if args.once:
            break
        time.sleep(15)
