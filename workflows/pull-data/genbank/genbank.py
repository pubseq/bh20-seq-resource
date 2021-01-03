# Genbank XML parser

from collections import namedtuple
import dateutil
from dateutil.parser import parse as dateparse
from dateutil.tz import gettz
import os
import re
import sys
import types
import xml.etree.ElementTree as ET

class GBError(Exception):
    pass

"""
Example of an output JSON:

{
  "id": "placeholder",
  "host": {
    "host_species": "http://purl.obolibrary.org/obo/NCBITaxon_9606"
  },
  "sample": {
    "sample_id": "MT890462.1",
    "source_database_accession": [
      "http://identifiers.org/insdc/MT890462.1#sequence"
    ],
    "collection_location": "http://www.wikidata.org/entity/Q649",
    "collection_date": "2020-04-17",
    "collecting_institution": "N.A.Kovtun Clinical Hospital 1 of Departament of President Affairs"
  },
  "virus": {
    "virus_strain": "SARS-CoV-2/human/RUS/20200417_10/2020",
    "virus_species": "http://purl.obolibrary.org/obo/NCBITaxon_2697049"
  },
  "technology": {
    "assembly_method": "http://purl.obolibrary.org/obo/GENEPIO_0001628",
    "alignment_protocol": "bowtie2 v. 2.3.4",
    "sample_sequencing_technology": [
      "http://purl.obolibrary.org/obo/OBI_0000759"
    ]
  },
  "submitter": {
    "authors": [
      "Blagodatskikh,K.A."
    ],
    "submitter_name": [
      "R&D"
    ],
    "submitter_address": "Pirogov Russian National Research Medical University, Ostrovityanova 1, Moscow 117997, Russia"
  }
}

Note: missing data should be None! Do not fill in other data by
'guessing'.

"""

def get_metadata(id, gbseq):
    host = types.SimpleNamespace()
    sample = types.SimpleNamespace()
    submitter = types.SimpleNamespace()
    warnings = []

    def warn(msg):
        print(f"WARNING: {msg}",file=sys.stderr)
        warnings.append(msg)

    host.host_species = "http://purl.obolibrary.org/obo/NCBITaxon_9606"
    sample.sample_id = id
    sample.database = "https://www.ncbi.nlm.nih.gov/genbank/"
    sample.source_database_accession = f"http://identifiers.org/insdc/{id}#sequence"
    # <GBQualifier>
    #   <GBQualifier_name>country</GBQualifier_name>
    #   <GBQualifier_value>USA: Cruise_Ship_1, California</GBQualifier_value>
    # </GBQualifier>
    sample.collection_location = "FIXME"

    submitter.authors = [n.text for n in gbseq.findall(".//GBAuthor")]
    # <GBReference_journal>Submitted (28-OCT-2020) MDU-PHL, The Peter
    #   Doherty Institute for Infection and Immunity, 792 Elizabeth
    #   Street, Melbourne, Vic 3000, Australia
    # </GBReference_journal>
    try:
        n = gbseq.find(".//GBReference_journal").text
        # print(n,file=sys.stderr)
        if n != 'Unpublished':
            institute,address = n.split(',',1)
            submitter.submitter_name = institute.split(') ')[1]
            submitter.submitter_address = address.strip()
    except AttributeError:
        pass
    except ValueError:
        submitter.additional_submitter_information = n
        pass

    # --- Dates
    n = gbseq.find("./GBSeq_create-date")
    creation_date = dateparse(n.text).date()
    n = gbseq.find("./GBSeq_update-date")
    update_date = dateparse(n.text).date()
    n = gbseq.find(".//GBQualifier/GBQualifier_name/[.='collection_date']/../GBQualifier_value")
    try:
        date = dateparse(n.text).date()
        sample.collection_date = str(date)
    except dateutil.parser._parser.ParserError as e:
        warn("No collection_date: ",str(e))
        sample.collection_date = None
    except AttributeError:
        warn("Missing collection_date")
        sample.collection_date = None

    info = {
        'id': 'placeholder',
        'update_date': str(update_date),
        'host': host,
        'sample': sample,
        #'virus': virus,
        #'technology': technology,
        'submitter': submitter,
        'warnings': warnings,
        }
    print(info)
    return True,info

def get_sequence(id, gbseq):
    seq = None
    count = 0
    for gbseq_sequence in gbseq.findall('./GBSeq_sequence'):
        count += 1
        if count > 1:
            raise GBError(f"Expected one sequence for {id}")
        seq = gbseq_sequence.text.upper()
        print(f"SEQ: size={len(seq)}",seq[0:30])
        if len(seq) < 20_000:
            raise GBError(f"Sequence too short")
        return seq
