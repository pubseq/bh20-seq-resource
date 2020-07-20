# Public API for PubSeq

import sys
import requests

from flask import Flask, request, redirect, send_file, send_from_directory, render_template, jsonify
from bh20simplewebuploader.main import app

@app.route('/api/version')
def version():
    return jsonify({ 'service': 'PubSeq', 'version': 0.10 })

@app.route('/api/ebi/sample-<id>.xml', methods=['GET'])
def ebi_sample(id):
    page = render_template('ebi-sample.xml',**locals())
    return page
