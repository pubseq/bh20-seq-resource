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
ARVADOS_API_TOKEN='2fbebpmbo3rw3x05ueu2i6nx70zhrsb1p22ycu3ry34m4x4462'
UPLOAD_PROJECT='lugli-j7d0g-n5clictpuvwk8aa'

def qa_stuff(metadata, sequence_p1, sequence_p2):
    try:
        log.debug("Checking metadata")
        if not qc_metadata(metadata.name):
            log.warning("Failed metadata qc")
            exit(1)
    except ValueError as e:
        log.debug(e)
        log.debug("Failed metadata qc")
        print(e)
        exit(1)

    target = []
    try:
        log.debug("Checking FASTA/FASTQ QC")
        target.append(qc_fasta(sequence_p1))
        if sequence_p2:
            target.append(qc_fasta(sequence_p2))
            target[0] = ("reads_1."+target[0][0][6:], target[0][1])
            target[1] = ("reads_2."+target[1][0][6:], target[0][1])
    except ValueError as e:
        log.debug(e)
        log.debug("Failed FASTA qc")
        print(e)
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
    parser.add_argument('sequence_p2', type=argparse.FileType('rb'), default=None, help='sequence FASTQ pair')
    parser.add_argument("--validate", action="store_true", help="Dry run, validate only")
    args = parser.parse_args()

    api = arvados.api(host=ARVADOS_API_HOST, token=ARVADOS_API_TOKEN, insecure=True)

    target = qa_stuff(args.metadata, args.sequence_p1, args.sequence_p2)
    seqlabel = target[0][1]

    if args.validate:
        print("Valid")
        exit(0)

    col = arvados.collection.Collection(api_client=api)

    upload_sequence(col, target[0], args.sequence_p1)
    if args.sequence_p2:
        upload_sequence(col, target[1], args.sequence_p2)

    print("Reading metadata")
    with col.open("metadata.yaml", "w") as f:
        r = args.metadata.read(65536)
        print(r[0:20])
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

    col.save_new(owner_uuid=UPLOAD_PROJECT, name="%s uploaded by %s from %s" %
                 (seqlabel, properties['upload_user'], properties['upload_ip']),
                 properties=properties, ensure_unique_name=True)

    print("Saved to %s" % col.manifest_locator())

    print("Done")

if __name__ == "__main__":
    main()
