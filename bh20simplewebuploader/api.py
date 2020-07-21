# Public API for PubSeq

import os
import requests
import sys

from flask import Flask, request, redirect, send_file, send_from_directory, render_template, jsonify
from bh20simplewebuploader.main import app, baseURL

@app.route('/api/version')
def version():
    return jsonify({ 'service': 'PubSeq', 'version': 0.10 })

@app.route('/api/ebi/sample-<id>.xml', methods=['GET'])
def ebi_sample(id):
    page = render_template('ebi-sample.xml',**locals())
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
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    # metadata = file.name(seq)+"/metadata.yaml"
    print(result)
    return jsonify([{
        'id': x['id']['value'],
        'fasta': x['seq']['value'],
        'collection': os.path.dirname(x['seq']['value']),
        'info': x['info']['value'],
    } for x in result])
