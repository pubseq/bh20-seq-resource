import argparse
import time
import arvados
import arvados.collection
import json
import urllib.request
import socket
import getpass
from .qc_metadata import qc_metadata

ARVADOS_API_HOST='lugli.arvadosapi.com'
ARVADOS_API_TOKEN='2fbebpmbo3rw3x05ueu2i6nx70zhrsb1p22ycu3ry34m4x4462'
UPLOAD_PROJECT='lugli-j7d0g-n5clictpuvwk8aa'

def main():
    parser = argparse.ArgumentParser(description='Upload SARS-CoV-19 sequences for analysis')
    parser.add_argument('sequence', type=argparse.FileType('r'), help='sequence FASTA')
    parser.add_argument('metadata', type=argparse.FileType('r'), help='sequence metadata json')
    args = parser.parse_args()

    api = arvados.api(host=ARVADOS_API_HOST, token=ARVADOS_API_TOKEN, insecure=True)

    if not qc_metadata(args.metadata.name):
        print("Failed metadata qc")
        exit(1)

    col = arvados.collection.Collection(api_client=api)

    if args.sequence.name.endswith("fasta") or args.sequence.name.endswith("fa"):
        target = "sequence.fasta"
    elif args.sequence.name.endswith("fastq") or args.sequence.name.endswith("fq"):
        target = "reads.fastq"

    with col.open(target, "w") as f:
        r = args.sequence.read(65536)
        print(r[0:20])
        while r:
            f.write(r)
            r = args.sequence.read(65536)

    print("Reading metadata")
    with col.open("metadata.yaml", "w") as f:
        r = args.metadata.read(65536)
        print(r[0:20])
        while r:
            f.write(r)
            r = args.metadata.read(65536)

    external_ip = urllib.request.urlopen('https://ident.me').read().decode('utf8')

    properties = {
        "upload_app": "bh20-seq-uploader",
        "upload_ip": external_ip,
        "upload_user": "%s@%s" % (getpass.getuser(), socket.gethostname())
    }

    col.save_new(owner_uuid=UPLOAD_PROJECT, name="Uploaded by %s from %s" %
                 (properties['upload_user'], properties['upload_ip']),
                 properties=properties, ensure_unique_name=True)

    print("Done")

if __name__ == "__main__":
    main()
