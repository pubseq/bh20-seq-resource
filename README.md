# Sequence uploader

This repository provides a sequence uploader for the

# Run

Run the uploader with a FASTA file and accompanying metadata:

    python3 bh20sequploader/main.py example/sequence.fasta example/metadata.json

# Add a workflow

get your SARS-CoV-2 sequences from GenBank in seqs.fa

```sh
minimap2 -cx asm20 -X seqs.fa seqs.fa >seqs.paf
seqwish -s seqs.fa -p seqs.paf -g seqs.gfa
odgi build -g seqs.gfa -s -o seqs.odgi
odgi viz -i seqs.odgi -o seqs.png -x 4000 -y 500 -R -P 5
```

from https://github.com/virtual-biohackathons/covid-19-bh20/wiki/Pangenome#pangenome-model-from-available-genomes

# Installation

This tool requires the arvados Python module which can be installed
using .deb or .rpm packages through
https://doc.arvados.org/v2.0/sdk/python/sdk-python.html. The actual
code lives [here](https://github.com/arvados/arvados/tree/master/sdk/python) and
suggests a local install using

    apt-get install libcurl4-openssl-dev libssl1.0-dev
    pip3 install --user arvados-python-client

Next update

    export PATH=$PATH:$HOME/.local/bin

## Install with GNU Guix

Set up a container:

    ~/opt/guix/bin/guix environment -C guix --ad-hoc python openssl python-pycurl nss-certs
    pip3 install --user arvados-python-client

Pip installed the following modules

    arvados-python-client-2.0.1 ciso8601-2.1.3 future-0.18.2 google-api-python-client-1.6.7 httplib2-0.17.1 oauth2client-4.1.3 pyasn1-0.4.8 pyasn1-modules-0.2.8 rsa-4.0 ruamel.yaml-0.15.77 six-1.14.0 uritemplate-3.0.1 ws4py-0.5.1
