# Public API for PubSeq

import os
import requests
import sys

from flask import Flask, request, redirect, send_file, send_from_directory, render_template, jsonify
from bh20simplewebuploader.main import app, sparqlURL

PUBSEQ="http://covid19.genenetwork.org"
ARVADOS="https://collections.lugli.arvadosapi.com/c="

# Helper functions

def fetch_sample_metadata(id):
    query = """
    PREFIX pubseq: <http://biohackathon.org/bh20-seq-schema#MainSchema/>
    PREFIX sio: <http://semanticscience.org/resource/>
    PREFIX edam: <http://edamontology.org/>
    PREFIX efo: <http://www.ebi.ac.uk/efo/>
    PREFIX evs: <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#>
    PREFIX obo: <http://purl.obolibrary.org/obo/>
    select distinct ?id ?seq ?date ?info ?specimen ?sequencer ?mapper
    {
      ?sample sio:SIO_000115 "%s" ;
              sio:SIO_000115 ?id ;
              evs:C25164 ?date .
      ?seq    pubseq:technology ?tech ;
              pubseq:sample ?sample .
      ?tech   efo:EFO_0002699 ?mapper ;
              obo:OBI_0600047 ?sequencer .
      optional { ?sample edam:data_2091 ?info } .
      optional { ?sample obo:OBI_0001479 ?specimen } .
    } limit 5
    """ % id
    payload = {'query': query, 'format': 'json'}
    r = requests.get(sparqlURL, params=payload)
    return r.json()['results']['bindings']

# Main API routes

@app.route('/api/version')
def version():
    return jsonify({ 'service': 'PubSeq', 'version': 0.10 })

@app.route('/api/sample/<id>.json')
def sample(id):
    """

API sample should return a record pointing to other resources,
notably: permalink, original metadata record and the fasta
data.

curl http://localhost:5067/api/sample/MT533203.1.json
[
  {
    "collection": "http://covid19.genenetwork.org/resource/lugli-4zz18-uovend31hdwa5ks",
    "date": "2020-04-27",
    "fasta": "https://collections.lugli.arvadosapi.com/c=lugli-4zz18-uovend31hdwa5ks/sequence.fasta",
    "id": "MT533203.1",
    "info": "http://identifiers.org/insdc/MT533203.1#sequence",
    "mapper": "minimap v. 2.17",
    "metadata": "https://collections.lugli.arvadosapi.com/c=lugli-4zz18-uovend31hdwa5ks/metadata.yaml",
    "permalink": "http://covid19.genenetwork.org/resource/MT533203.1",
    "sequencer": "http://www.ebi.ac.uk/efo/EFO_0008632",
    "specimen": "http://purl.obolibrary.org/obo/NCIT_C155831"
  }
]


"""
    # metadata = file.name(seq)+"/metadata.yaml"
    meta = fetch_sample_metadata(id)
    print(meta)
    # http://collections.lugli.arvadosapi.com/c=lugli-4zz18-uovend31hdwa5ks/metadata.yaml
    return jsonify([{
        'id': x['id']['value'],
        'collection': x['seq']['value'],
        'permalink': PUBSEQ+'/resource/'+x['id']['value'],
        'fasta': ARVADOS+os.path.basename(x['seq']['value'])+'/sequence.fasta',
        'metadata': ARVADOS+os.path.basename(x['seq']['value'])+'/metadata.yaml',
        'date': x['date']['value'],
        'info': x['info']['value'],
        'specimen': x['specimen']['value'],
        'sequencer': x['sequencer']['value'],
        'mapper': x['mapper']['value'],
    } for x in meta])

@app.route('/api/ebi/sample-<id>.xml', methods=['GET'])
def ebi_sample(id):
    meta = fetch_sample_metadata(id)[0]
    page = render_template('ebi-sample.xml',sampleid=id,sequencer=meta['sequencer']['value'],date=meta['date']['value'],specimen=meta['specimen']['value'])
    return page

@app.route('/api/search', methods=['GET'])
def search():
    """
    Execute a 'global search'
    """
    s = request.args.get('s')
    if s == "":
        s = "MT326090.1"
    query = """
    PREFIX pubseq: <http://biohackathon.org/bh20-seq-schema#MainSchema/>
    PREFIX sio: <http://semanticscience.org/resource/>
    PREFIX edam: <http://edamontology.org/>
    select distinct ?id ?seq ?info
    {
    ?sample sio:SIO_000115 "%s" .
    ?sample sio:SIO_000115 ?id .
    ?seq pubseq:sample ?sample .
    ?sample edam:data_2091 ?info .
    } limit 100
    """ % s
    payload = {'query': query, 'format': 'json'}
    r = requests.get(sparqlURL, params=payload)
    result = r.json()['results']['bindings']
    # metadata = file.name(seq)+"/metadata.yaml"
    print(result)
    return jsonify([{
        'id': x['id']['value'],
        'fasta': x['seq']['value'],
        'collection': os.path.dirname(x['seq']['value']),
        'info': x['info']['value'],
    } for x in result])
