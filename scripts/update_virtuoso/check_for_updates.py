#! /usr/bin/env python3
#
# Check for updates on Arvados, pull the TTL and
# push into Virtuoso
#
# You can run this in a Guix container with
#
#  ~/opt/guix/bin/guix environment -C guix --ad-hoc python python-requests curl --network -- python3 ./scripts/update_virtuoso/check_for_updates.py cache.txt dba dba

import requests
import time
import sys

assert(len(sys.argv)==4)
fn = sys.argv[1]
user = sys.argv[2]
pwd = sys.argv[3]


def upload(fn):
    # Upload into Virtuoso using CURL
    # cmd = "curl -X PUT --digest -u dba:dba -H Content-Type:text/turtle -T metadata.ttl -G http://localhost:8890/sparql-graph-crud-auth --data-urlencode graph=http://covid-19.genenetwork.org/graph".split(" ")
    # print("DELETE "+fn)
    # cmd = ("curl --digest --user dba:%s --verbose --url -G http://sparql.genenetwork.org/sparql-graph-crud-auth --data-urlencode graph=http://covid-19.genenetwork.org/graph -X DELETE" % pwd).split(" ")

    print("UPLOAD "+fn)
    cmd = ("curl -X PUT --digest -u dba:%s -H Content-Type:text/turtle -T %s -G http://sparql.genenetwork.org/sparql-graph-crud-auth --data-urlencode graph=http://covid-19.genenetwork.org/graph" % (pwd, fn) ).split(" ")
    print(cmd)
    p = subprocess.Popen(cmd)
    output = p.communicate()[0]
    print(output)
    assert(p.returncode == 0)

url = 'https://download.lugli.arvadosapi.com/c=lugli-4zz18-z513nlpqm03hpca/_/mergedmetadata.ttl'
# --- Fetch headers from TTL file on Arvados
r = requests.head(url)
print(r.headers)

print(r.headers['Last-Modified'])

# --- Convert/validate time stamp
# ValueError: time data 'Tue, 21 Apr 2020 23:47:43 GMT' does not match format '%a %b %d %H:%M:%S %Y'
last_modified_str = r.headers['Last-Modified']
t_stamp = time.strptime(last_modified_str,"%a, %d %b %Y %H:%M:%S %Z" )
print(t_stamp)

# OK, it works, now check last stored value
import os.path
stamp = None
if os.path.isfile(fn):
    file = open(fn,"r")
    stamp = file.read()
    file.close

import subprocess
if stamp != last_modified_str:
    print("Fetch metadata TTL")
    r = requests.get(url)
    assert(r.status_code == 200)
    with open("metadata.ttl", "w") as f:
        f.write(r.text)
        f.close
    upload("metadata.ttl")
    upload("semantic_enrichment/labels.ttl")
    upload("semantic_enrichment/countries.ttl")

    with open(fn,"w") as f:
        f.write(last_modified_str)
