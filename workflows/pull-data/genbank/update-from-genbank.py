#!/usr/bin/env python3
#
# - bulk download genbank data and matadata, preparing the FASTA and
#   the YAML files
#
# See .guix-run

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ids-to-ignore', type=str, help='file with ids to ignore in all steps, 1 id per line', required=False)
parser.add_argument('--ids-to-consider', type=str, help='file with ids to consider in all steps, 1 id per line', required=False)
parser.add_argument('--skip-request', action='store_true', help='skip metadata and sequence request', required=False)
parser.add_argument('--only-missing-ids', action='store_true', help='download only missing ids not already downloaded', required=False)
parser.add_argument('--dict-ontology', type=str, help='where is the ontology',
                    default='../dict_ontology_standardization/', required=False)
args = parser.parse_args()

from Bio import Entrez
Entrez.email = 'another_email@gmail.com'

import xml.etree.ElementTree as ET
import json
import os
import requests

from datetime import date, datetime
from dateutil.parser import parse

import sys
sys.path.append('../')
from utils import is_integer, chunks, check_and_get_ontology_dictionaries


num_ids_for_request = 100
min_acceptable_collection_date = datetime(2019, 12, 1)

dir_metadata = 'metadata_from_nuccore'
dir_fasta_and_yaml = 'fasta_and_yaml'
dir_dict_ontology_standardization = args.dict_ontology

today_date = date.today().strftime("%Y.%m.%d")
path_ncbi_virus_accession = 'sequences.{}.acc'.format(today_date)


field_to_term_to_uri_dict = check_and_get_ontology_dictionaries(dir_dict_ontology_standardization)


if os.path.exists(dir_metadata):
    print("The directory '{}' already exists.".format(dir_metadata))

    if not args.skip_request:
        print("\tTo start the request, delete the directory '{}' or specify --skip-request.".format(dir_metadata))
        sys.exit(-1)


accession_to_ignore_set = set()

if args.ids_to_ignore:
    if not os.path.exists(args.ids_to_ignore):
        print("\tThe '{}' file doesn't exist.".format(args.ids_to_ignore))
        sys.exit(-1)

    with open(args.ids_to_ignore) as f:
        accession_to_ignore_set.update(set([x.split('.')[0] for x in f.read().strip('\n').split('\n')]))

    print('There are {} accessions to ignore.'.format(len(accession_to_ignore_set)))


# ----------------------------------------------------------------------
"""
With --only-missing-ids only download accessions that we do not yet have!
"""
accession_already_downloaded_set = set()

if os.path.exists(dir_fasta_and_yaml):
    """
    If the fasta_and_yaml directory exists and --only-missing-ids was set
    we make a list of all downloaded accessions:
    """
    print("The directory '{}' already exists.".format(dir_fasta_and_yaml))
    if not args.only_missing_ids:
        print("To start the download, delete the directory '{}' or specify --only-missing-ids.".format(dir_fasta_and_yaml))
        sys.exit(-1)

    """
    Fetch all YAML filenames and load `accession_already_downloaded_set`
    """
    accession_already_downloaded_set = set([x.split('.yaml')[0].split('.')[0] for x in os.listdir(dir_fasta_and_yaml) if x.endswith('.yaml')])
    print('There are {} accessions already downloaded.'.format(len(accession_already_downloaded_set)))

accession_to_ignore_set.update(accession_already_downloaded_set)

# ----------------------------------------------------------------------
"""
Check for --ids-to-consider
"""
accession_to_consider_set = set()

if args.ids_to_consider:
    if not os.path.exists(args.ids_to_consider):
        print("\tThe '{}' file doesn't exist.".format(args.ids_to_consider))
        sys.exit(-1)

    with open(args.ids_to_consider) as f:
        accession_to_consider_set.update(set([x.split('.')[0] for x in f.read().strip('\n').split('\n')]))

    if len(accession_to_consider_set) > 0:
        print('There are {} accessions to consider.'.format(len(accession_to_consider_set)))

# ----------------------------------------------------------------------
"""
Download section for genbank XML
"""

if not os.path.exists(dir_metadata):
    # Take all the ids
    id_set = set()

    # Try to search several strings
    term_list = ['SARS-CoV-2', 'SARS-CoV2', 'SARS CoV2', 'SARSCoV2', 'txid2697049[Organism]']
    for term in term_list:
        tmp_list = Entrez.read(
            Entrez.esearch(db='nuccore', term=term, idtype='acc', retmax='10000')
        )['IdList']

        # Remove mRNAs, ncRNAs, Proteins, and predicted models (more information here: https://en.wikipedia.org/wiki/RefSeq)
        # Remove the version in the id
        new_ids_set = set([x.split('.')[0] for x in tmp_list if x[:2] not in ['NM', 'NR', 'NP', 'XM', 'XR', 'XP', 'WP']])

        if len(accession_to_consider_set) > 0:
            new_ids_set = new_ids_set.intersection(accession_to_consider_set)

        new_ids = len(new_ids_set.difference(id_set))
        id_set.update(new_ids_set)

        print('Term:', term, '-->', new_ids, 'new IDs from', len(tmp_list), '---> Total unique IDs:', len(id_set))

    if not os.path.exists(path_ncbi_virus_accession):
        r = requests.get('https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?q=*:*&fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(2697049)&cmd=download&sort=SourceDB_s%20desc,CreateDate_dt%20desc,id%20asc&dlfmt=acc&fl=id')
        with open(path_ncbi_virus_accession, 'w') as fw:
            fw.write(r.text)

    with open(path_ncbi_virus_accession) as f:
        tmp_list = [line.strip('\n') for line in f]

    new_ids_set = set(tmp_list)
    if len(accession_to_consider_set) > 0:
        new_ids_set = new_ids_set.intersection(accession_to_consider_set)

    new_ids = len(new_ids_set.difference(id_set))
    id_set.update(new_ids_set)

    print('DB: NCBI Virus', today_date, '-->', new_ids, 'new IDs from', len(tmp_list), '---> Total unique IDs:', len(id_set))

    id_set = id_set.difference(accession_to_ignore_set)
    print('There are {} missing IDs to download.'.format(len(id_set)))

    os.makedirs(dir_metadata)
    for i, id_x_list in enumerate(chunks(list(id_set), num_ids_for_request)):
        path_metadata_xxx_xml = os.path.join(dir_metadata, 'metadata_{}.xml'.format(i))
        print('Requesting {} ids --> {}'.format(len(id_x_list), path_metadata_xxx_xml))

        with open(path_metadata_xxx_xml, 'w') as fw:
            fw.write(
                Entrez.efetch(db='nuccore', id=id_x_list, retmode='xml').read()
            )
