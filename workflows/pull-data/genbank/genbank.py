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
    "collecting_institution": "N.A.Kovtun Clinical Hospital 1 of Departament of President Affairs",
    "specimen_source": ["http://purl.obolibrary.org/obo/NCIT_C155831"]
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


"""

def get_metadata(id, gbseq):
    """This is a minimal data parser from genbank XML records. Inference
    on, for example geo location, is not allowed in this function and
    happens downstream.

    That is to keep the parsing simple.

    Important: missing data should be missing or None! Do not fill in
    data by 'guessing'.

    When data is malformed a warning should be logged and added to the
    warning list.

    """
    host = types.SimpleNamespace()
    sample = types.SimpleNamespace()
    submitter = types.SimpleNamespace()
    technology = types.SimpleNamespace()
    virus = types.SimpleNamespace()
    warnings = []

    def warn(msg):
        print(f"WARNING: {msg}",file=sys.stderr)
        warnings.append(msg)

    def fetch(msg, xpath):
        try:
            n = gbseq.find(xpath).text
            return n
        except AttributeError:
            warn("Missing "+msg)

    host.host_species = "http://purl.obolibrary.org/obo/NCBITaxon_9606"
    sample.sample_id = id
    sample.database = "https://www.ncbi.nlm.nih.gov/genbank/"
    sample.source_database_accession = f"http://identifiers.org/insdc/{id}#sequence"
    #   <GBQualifier_value>USA: Cruise_Ship_1, California</GBQualifier_value>
    n = fetch("host_species", ".//GBQualifier/GBQualifier_name/[.='country']/../GBQualifier_value")
    if n: sample.collection_location = n
    else: warn("Missing collection_location")

    submitter.authors = [n.text for n in gbseq.findall(".//GBAuthor")]
    if not len(submitter.authors): warn("Missing authors")

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

    try:
        n = gbseq.find("./GBSeq_comment").text
    except AttributeError:
        pass
    if 'Assembly-Data' in n:
        # print(n,file=sys.stderr)
        # the following is wrong (de novo by default)
        # technology.assembly_method = 'http://purl.obolibrary.org/obo/GENEPIO_0001628'
        p = re.compile(r'.*Assembly Method :: ([^;]+).*')
        m = p.match(n)
        if m: technology.alignment_protocol = m.group(1)
        p = re.compile(r'.*Coverage :: ([^;]+).*')
        m = p.match(n)
        if m: technology.sequencing_coverage = m.group(1)
        p = re.compile(r'.*Sequencing Technology :: ([^;]+).*')
        m = p.match(n)
        if m: technology.sample_sequencing_technology = m.group(1).strip()
        else: warn("Missing sample_sequencing_technology")

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
    except AttributeError:
        warn("Missing collection_date")

    # --- Host info
    # - Homo sapiens
    # - Homo sapiens; female
    # - Homo sapiens; female 63
    # - Homo sapiens; female; age 40
    # - Homo sapiens; gender: F; age: 61
    # - Homo sapiens; gender: M; age: 68
    # - Homo sapiens; hospitalized patient
    # - Homo sapiens; male
    # - Homo sapiens; male; 63
    # - Homo sapiens; male; age 29
    # - Homo sapiens; symptomatic
    n = fetch("host_species", ".//GBQualifier/GBQualifier_name/[.='host']/../GBQualifier_value")
    if n:
        list = n.split('; ')
        species = list[0]
        host.host_species = species
        if species != "Homo sapiens":
            warn(f"Species not understood: {species}")
        if len(list)>1:
            sex = list[1]
            if 'male' in sex or 'gender: M' in sex: host.host_sex = 'male'
            if 'female' in sex or 'gender: F' in sex: host.host_sex = 'female'
        if len(list)>2:
            age = list[2]
            p = re.compile(r'[^\d]+(\d+)')
            m = p.match(n)
            print(m.group(1))
            if m:
                host.host_age = int(m.group(1))
                host.host_age_unit = 'http://purl.obolibrary.org/obo/UO_0000036'
        # sys.exit(1)
    n = fetch("virus_strain", ".//GBQualifier/GBQualifier_name/[.='isolate']/../GBQualifier_value")
    if n: virus.virus_strain = n
    n = fetch("virus_species", ".//GBQualifier/GBQualifier_name/[.='db_xref']/../GBQualifier_value")
    if n: virus.virus_species = "http://purl.obolibrary.org/obo/NCBITaxon_"+n.split('taxon:')[1]
    n = fetch("specimen_source", ".//GBQualifier/GBQualifier_name/[.='isolation_source']/../GBQualifier_value")
    if n: sample.specimen_source = n

    info = {
        'id': 'placeholder',
        'update_date': str(update_date),
        'host': host.__dict__,
        'sample': sample.__dict__,
        'virus': virus.__dict__,
        'technology': technology.__dict__,
        'submitter': submitter.__dict__,
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
