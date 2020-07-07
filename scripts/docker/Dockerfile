FROM debian:10

RUN apt-get update && \
    apt-get -yq --no-install-recommends -o Acquire::Retries=6 install \
    python3 python3-pip python3-setuptools python3-dev python-pycurl \
    minimap2 python3-biopython libcurl4-openssl-dev build-essential \
    libssl-dev libmagic-dev python3-magic && \
    apt-get clean

RUN pip3 install bh20-seq-uploader py-dateutil
