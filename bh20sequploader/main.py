import argparse
import time
import arvados
import arvados.collection
import json
import logging
import magic
from pathlib import Path
import urllib.request
import socket
import getpass
import sys
sys.path.insert(0,'.')
from bh20sequploader.qc_metadata import qc_metadata
from bh20sequploader.qc_fasta import qc_fasta

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__ )
log.debug("Entering sequence uploader")

ARVADOS_API_HOST='lugli.arvadosapi.com'
UPLOADER_API_TOKEN='2fbebpmbo3rw3x05ueu2i6nx70zhrsb1p22ycu3ry34m4x4462'
ANONYMOUS_API_TOKEN='5o42qdxpxp5cj15jqjf7vnxx5xduhm4ret703suuoa3ivfglfh'
UPLOAD_PROJECT='lugli-j7d0g-n5clictpuvwk8aa'
VALIDATED_PROJECT='lugli-j7d0g-5ct8p1i1wrgyjvp'

def qc_stuff(metadata, sequence_p1, sequence_p2, do_qc=True):
    failed = False
    try:
        log.debug("Checking metadata" if do_qc else "Skipping metadata check")
        if do_qc and not qc_metadata(metadata.name):
            log.warning("Failed metadata QC")
            failed = True
    except Exception as e:
        log.exception("Failed metadata QC")
        failed = True

    target = []
    try:
        log.debug("FASTA/FASTQ QC" if do_qc else "Limited FASTA/FASTQ QC")
        target.append(qc_fasta(sequence_p1, check_with_mimimap2=do_qc))
        if sequence_p2:
            target.append(qc_fasta(sequence_p2))
            target[0] = ("reads_1."+target[0][0][6:], target[0][1])
            target[1] = ("reads_2."+target[1][0][6:], target[0][1])
    except Exception as e:
        log.exception("Failed sequence QC")
        failed = True

    if failed:
        exit(1)

    return target

def upload_sequence(col, target, sequence):
    with col.open(target[0], "wb") as f:
        r = sequence.read(65536)
        while r:
            f.write(r)
            r = sequence.read(65536)


def main():
    parser = argparse.ArgumentParser(description='Upload SARS-CoV-19 sequences for analysis')
    parser.add_argument('metadata', type=argparse.FileType('r'), help='sequence metadata json')
    parser.add_argument('sequence_p1', type=argparse.FileType('rb'), help='sequence FASTA/FASTQ')
    parser.add_argument('sequence_p2', type=argparse.FileType('rb'), default=None, nargs='?', help='sequence FASTQ pair')
    parser.add_argument("--validate", action="store_true", help="Dry run, validate only")
    parser.add_argument("--skip-qc", action="store_true", help="Skip local qc check")
    parser.add_argument("--trusted", action="store_true", help="Trust local validation and add directly to validated project")
    args = parser.parse_args()

    if args.trusted:
        # Use credentials from environment
        api = arvados.api()
    else:
        api = arvados.api(host=ARVADOS_API_HOST, token=UPLOADER_API_TOKEN, insecure=True)

    target = qc_stuff(args.metadata, args.sequence_p1, args.sequence_p2, not args.skip_qc)
    seqlabel = target[0][1]

    if args.validate:
        log.info("Valid")
        exit(0)

    col = arvados.collection.Collection(api_client=api)

    upload_sequence(col, target[0], args.sequence_p1)
    if args.sequence_p2:
        upload_sequence(col, target[1], args.sequence_p2)

    log.info("Reading metadata")
    with col.open("metadata.yaml", "w") as f:
        r = args.metadata.read(65536)
        log.info(r[0:20])
        while r:
            f.write(r)
            r = args.metadata.read(65536)

    external_ip = urllib.request.urlopen('https://ident.me').read().decode('utf8')

    try:
        username = getpass.getuser()
    except KeyError:
        username = "unknown"

    properties = {
        "sequence_label": seqlabel,
        "upload_app": "bh20-seq-uploader",
        "upload_ip": external_ip,
        "upload_user": "%s@%s" % (username, socket.gethostname())
    }

    api2 = arvados.api(host=ARVADOS_API_HOST, token=ANONYMOUS_API_TOKEN, insecure=True)
    dup = api2.collections().list(filters=[["owner_uuid", "in", [VALIDATED_PROJECT, UPLOAD_PROJECT]],
                                           ["portable_data_hash", "=", col.portable_data_hash()]]).execute()
    if dup["items"]:
        # This exact collection has been uploaded before.
        log.error("Duplicate of %s" % ([d["uuid"] for d in dup["items"]]))
        exit(1)

    if args.trusted:
        properties["status"] = "validated"
        owner_uuid = VALIDATED_PROJECT
    else:
        owner_uuid = UPLOAD_PROJECT

    col.save_new(owner_uuid=owner_uuid, name="%s uploaded by %s from %s" %
                 (seqlabel, properties['upload_user'], properties['upload_ip']),
                 properties=properties, ensure_unique_name=True)

    log.info("Saved to %s" % col.manifest_locator())
    log.info("Done")
    exit(0)

if __name__ == "__main__":
    main()
