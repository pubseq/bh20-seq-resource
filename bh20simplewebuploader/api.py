# Public API for PubSeq

import os
import requests
import sys

from flask import Flask, request, redirect, send_file, send_from_directory, render_template, jsonify
from bh20simplewebuploader.main import app, sparqlURL

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
    # metadata = file.name(seq)+"/metadata.yaml"
    meta = fetch_sample_metadata(id)
    print(meta)
    return jsonify([{
        'id': x['id']['value'],
        'fasta': x['seq']['value'],
        'collection': os.path.dirname(x['seq']['value']),
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
