# Public API for PubSeq

import os
import requests
import sys
import types

from flask import Flask, request, redirect, send_file, send_from_directory, render_template, jsonify
from bh20simplewebuploader.main import app, sparqlURL

PUBSEQ="http://covid19.genenetwork.org"
ARVADOS="https://collections.lugli.arvadosapi.com/c="

# Helper functions

def fetch_sample(id, query=None):
    default_query = """

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
      optional { ?tech   efo:EFO_0002699 ?mapper } .
      optional { ?tech   obo:OBI_0600047 ?sequencer . }
      optional { ?sample edam:data_2091 ?info } .
      optional { ?sample obo:OBI_0001479 ?specimen } .
    } limit 5

    """ % id
    if not query: query = default_query
    print(query)
    payload = {'query': query, 'format': 'json'}
    r = requests.get(sparqlURL, params=payload)
    res = r.json()
    print(res)
    return res['results']['bindings'],res['head']['vars']

def fetch_one_sample(id, query=None):
    """Get the top sample and return a SimpleNamespace"""

    result,varlist = fetch_sample(id,query)
    h = {}
    row = result[0]
    for key in varlist:
        if key in row:
            h[key] = row[key]['value']
    print(h)
    h['arv_id'] = os.path.basename(h['seq'])
    return types.SimpleNamespace(**h)

def fetch_one_record(id):
    m = fetch_one_sample(id)
    arv_id = m.arv_id
    rec = { "id": id,
            'arv_id': arv_id,
            "permalink": PUBSEQ+'/resource/'+id,
            "collection": m.seq,
            'collection_date': m.date,
            'fasta': ARVADOS+arv_id+'/sequence.fasta',
            'metadata': ARVADOS+arv_id+'/metadata.yaml',
    }
    h = m.__dict__ # for optional items
    if 'mapper' in h: rec['mapper'] = m.mapper
    if 'sequencer' in h: rec['sequencer']= m.sequencer
    return rec

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
{
  "id": "MT533203.1",
  "permalink": "http://covid19.genenetwork.org/resource/MT533203.1",
  "collection": "http://covid19.genenetwork.org/resource/lugli-4zz18-uovend31hdwa5ks",
  "collection_date": "2020-04-27",
  "fasta": "https://collections.lugli.arvadosapi.com/c=lugli-4zz18-uovend31hdwa5ks/sequence.fasta",
  "metadata": "https://collections.lugli.arvadosapi.com/c=lugli-4zz18-uovend31hdwa5ks/metadata.yaml",
  "mapper": "minimap v. 2.17",
  "sequencer": "http://www.ebi.ac.uk/efo/EFO_0008632"
}

"""

    return jsonify([fetch_one_record(id)])

@app.route('/api/ebi/sample-<id>.xml', methods=['GET'])
def ebi_sample(id):
    meta,varlist = fetch_sample(id)[0]
    page = render_template('ebi-sample.xml',sampleid=id,sequencer=meta['sequencer']['value'],date=meta['date']['value'],specimen=meta['specimen']['value'])
    return page

@app.route('/api/search', methods=['GET'])
def search():
    """
    Execute a 'global search'. Currently just duplicates fetch one
    sample. Should be more flexible FIXME.
    """
    s = request.args.get('s')
    if s == "": s = "MT326090.1"
    return jsonify([fetch_one_record(s)])
