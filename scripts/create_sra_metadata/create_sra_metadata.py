#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--ids-to-ignore', type=str, help='file with ids to ignore in all steps, 1 id per line', required=False)
parser.add_argument('--dict-ontology', type=str, help='where is the ontology',
                    default='../dict_ontology_standardization/', required=False)

args = parser.parse_args()

import os
from dateutil.parser import parse
import xml.etree.ElementTree as ET
import json
import gzip
import sys

dir_yaml = 'yaml'

date = '2020.07.09'

# Query on SRA: 'txid2697049[Organism]' (https://www.ncbi.nlm.nih.gov/sra/?term=txid2697049%5BOrganism%5D)
# Query on SRA: 'txid2697049[Organism:noexp] NOT 0[Mbases ' (https://www.ncbi.nlm.nih.gov/sra/?term=txid2697049%5BOrganism:noexp%5D%20NOT%200[Mbases)
#         -> Send to -> File -> Full XML -> Create File
path_sra_metadata_xml = 'SraExperimentPackage.{}.xml.gz'.format(date)

dir_dict_ontology_standardization = args.dict_ontology
path_sra_study_accessions_txt = 'SRAStudyAccessions.{}.txt'.format(date)


def is_integer(string_to_check):
    try:
        int(string_to_check)
        return True
    except ValueError:
        return False


accession_to_ignore_set = set()

if args.ids_to_ignore:
    if not os.path.exists(args.ids_to_ignore):
        print("\tThe '{}' file doesn't exist.".format(args.ids_to_ignore))
        sys.exit(-1)

    with open(args.ids_to_ignore) as f:
        accession_to_ignore_set.update(set([x.split('.')[0] for x in f.read().strip('\n').split('\n')]))
        print('There are {} accessions to ignore.'.format(len(accession_to_ignore_set)))


term_to_uri_dict = {}

for path_dict_xxx_csv in [os.path.join(dir_dict_ontology_standardization, name_xxx_csv) for name_xxx_csv in os.listdir(dir_dict_ontology_standardization) if name_xxx_csv.endswith('.csv')]:
    print('Read {}'.format(path_dict_xxx_csv))

    with open(path_dict_xxx_csv) as f:
        for line in f:
            if len(line.split(',')) > 2:
                term, uri = line.strip('\n').split('",')
            else:
                term, uri = line.strip('\n').split(',')

            term = term.strip('"')

            if term in term_to_uri_dict:
                print('Warning: in the dictionaries there are more entries for the same term ({}).'.format(term))
                continue

            term_to_uri_dict[term] = uri


if not os.path.exists(dir_yaml):
    os.makedirs(dir_yaml)


sra_metadata_xml_file = gzip.open(path_sra_metadata_xml, 'r')
tree = ET.parse(sra_metadata_xml_file)
sra_metadata_xml_file.close()

EXPERIMENT_PACKAGE_SET = tree.getroot()

missing_value_list = []
not_created_accession_dict = {}

run_accession_set = set()
run_accession_to_downloadble_file_url_dict = {}

num_yaml_created = 0

for i, EXPERIMENT_PACKAGE in enumerate(EXPERIMENT_PACKAGE_SET):
    #print(i, EXPERIMENT_PACKAGE)

    # A general default-empty yaml could be read from the definitive one
    info_for_yaml_dict = {
        'id': 'placeholder',
        'host': {},
        'sample': {},
        'virus': {},
        'technology': {},
        'submitter': {}
    }

    RUN_SET = EXPERIMENT_PACKAGE.find('RUN_SET')
    RUN = RUN_SET.find('RUN')
    accession = RUN.attrib['accession']
    run_accession_set.add(accession)
    #print(accession)

    info_for_yaml_dict['sample']['sample_id'] = accession

    #SRAFiles = RUN.find('SRAFiles')
    #if SRAFiles is not None:
    #    url = SRAFiles.find('SRAFile').attrib['url']
    #    if 'sra-download.ncbi.nlm.nih.gov' in url:
    #        run_accession_to_downloadble_file_url_dict[accession] = url


    SAMPLE = EXPERIMENT_PACKAGE.find('SAMPLE')
    SAMPLE_ATTRIBUTE_list = SAMPLE.iter('SAMPLE_ATTRIBUTE')

    for SAMPLE_ATTRIBUTE in SAMPLE_ATTRIBUTE_list:
        VALUE = SAMPLE_ATTRIBUTE.find('VALUE')
        if VALUE is not None:
            TAG_text = SAMPLE_ATTRIBUTE.find('TAG').text
            VALUE_text = VALUE.text

            if TAG_text in ['host', 'host scientific name']:
                if VALUE_text.lower() in ['homo sapien', 'homosapiens']:
                    VALUE_text = 'Homo sapiens'

                if VALUE_text in term_to_uri_dict:
                    info_for_yaml_dict['host']['host_species'] = term_to_uri_dict[VALUE_text]
                else:
                    missing_value_list.append('\t'.join([accession, 'host_species', VALUE_text]))
            elif TAG_text in ['host_health_status', 'host health state']:
                if VALUE_text in term_to_uri_dict:
                    info_for_yaml_dict['host']['host_health_status'] = term_to_uri_dict[VALUE_text]
                elif VALUE_text.strip("'") not in ['missing', 'not collected', 'not provided']:
                    missing_value_list.append('\t'.join([accession, 'host_health_status', VALUE_text]))
            elif TAG_text in ['strain', 'isolate']:
                if VALUE_text.lower() not in ['not applicable', 'missing', 'na', 'unknown', 'not provided']:
                    value_to_insert = VALUE_text

                    if value_to_insert.lower() in ['homo sapien', 'homosapiens']:
                        value_to_insert = 'Homo sapiens'

                    if value_to_insert in term_to_uri_dict:
                        value_to_insert = term_to_uri_dict[value_to_insert]

                    if 'virus_strain' not in info_for_yaml_dict:
                        info_for_yaml_dict['virus']['virus_strain'] = value_to_insert
                    else:
                        info_for_yaml_dict['virus']['virus_strain'] += '; ' + value_to_insert
            elif TAG_text in ['isolation_source', 'isolation source host-associated']:
                if VALUE_text in term_to_uri_dict:
                    info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict[VALUE_text]]
                else:
                    if VALUE_text.lower() in ['np/op', 'np/op swab', 'np/np swab', 'nasopharyngeal and oropharyngeal swab', 'nasopharyngeal/oropharyngeal swab', 'combined nasopharyngeal and oropharyngeal swab', 'naso and/or oropharyngeal swab']:
                        info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasopharyngeal swab'], term_to_uri_dict['oropharyngeal swab']]
                    elif VALUE_text.lower() in ['nasopharyngeal swab/throat swab', 'nasopharyngeal/throat swab', 'nasopharyngeal swab and throat swab', 'nasal swab and throat swab', 'nasopharyngeal aspirate/throat swab', 'Nasopharyngeal/Throat']:
                        info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasopharyngeal swab'], term_to_uri_dict['throat swab']]
                    elif VALUE_text.lower() in ['nasopharyngeal aspirate & throat swab', 'nasopharyngeal aspirate and throat swab']:
                        info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasopharyngeal aspirate'], term_to_uri_dict['throat swab']]
                    elif VALUE_text.lower() in ['nasal swab and throat swab']:
                        info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasal swab'], term_to_uri_dict['throat swab']]
                    elif VALUE_text.lower() in ['nasal-swab and oro-pharyngeal swab']:
                        info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasal swab'], term_to_uri_dict['oropharyngeal swab']]
                    elif VALUE_text.strip("'") not in ['missing', 'not collected', 'unknown', 'not provided', 'not applicable', 'N/A']:
                        missing_value_list.append('\t'.join([accession, 'specimen_source', VALUE_text]))
            elif TAG_text in ['host_sex', 'host sex']:
                if VALUE_text.lower() not in ['missing', 'not provided']:
                    if VALUE_text in ['male', 'female']:
                        info_for_yaml_dict['host']['host_sex'] = "http://purl.obolibrary.org/obo/PATO_0000384" if VALUE_text == 'male' else "http://purl.obolibrary.org/obo/PATO_0000383"
                    else:
                        missing_value_list.append('\t'.join([accession, 'host_sex', VALUE_text]))
            elif TAG_text in ['host_age', 'host age']:
                if is_integer(VALUE_text):
                    info_for_yaml_dict['host']['host_age'] = VALUE_text
                    info_for_yaml_dict['host']['host_age_unit'] = 'http://purl.obolibrary.org/obo/UO_0000036'
            elif TAG_text == 'collected_by':
                if VALUE_text.lower() not in ['not available', 'missing']:
                    name = VALUE_text in ['Dr. Susie Bartlett', 'Ahmed Babiker', 'Aisi Fu', 'Brandi Williamson', 'George Taiaroa', 'Natacha Ogando', 'Tim Dalebout', 'ykut Ozdarendeli']

                    info_for_yaml_dict['sample']['collector_name' if name else 'collecting_institution'] = VALUE_text
            elif TAG_text == 'collecting institution':
                if VALUE_text.lower() not in ['not provided', 'na']:
                    info_for_yaml_dict['sample']['collecting_institution'] = VALUE_text
            elif TAG_text in ['collection_date', 'collection date']:
                if VALUE_text.lower() not in ['not applicable', 'missing', 'na']:
                    date_to_write = VALUE_text
                    date_is_estimated = True

                    VALUE_text_list = VALUE_text.split('-')
                    if len(VALUE_text_list) == 3:
                        date_is_estimated = False

                        if VALUE_text_list[1].isalpha():
                            date_to_write = parse(VALUE_text).strftime('%Y-%m-%d')
                    elif len(VALUE_text_list) == 2:
                        date_to_write = VALUE_text + '-15'
                    else:
                        if int(VALUE_text) < 2020:
                            date_to_write = "{}-12-15".format(VALUE_text)
                        else:
                            date_to_write = "{}-01-15".format(VALUE_text)

                    info_for_yaml_dict['sample']['collection_date'] = date_to_write

                    if date_is_estimated:
                        if 'additional_collection_information' in info_for_yaml_dict['sample']:
                            info_for_yaml_dict['sample']['additional_collection_information'] += "; The 'collection_date' is estimated (the original date was: {})".format(VALUE_text)
                        else:
                            info_for_yaml_dict['sample']['additional_collection_information'] = "The 'collection_date' is estimated (the original date was: {})".format(VALUE_text)
            elif TAG_text in ['geo_loc_name', 'geographic location (country and/or sea)', 'geographic location (region and locality)']:
                if ': ' in VALUE_text:
                    VALUE_text = VALUE_text.replace(': ', ':')

                if VALUE_text in term_to_uri_dict:
                    info_for_yaml_dict['sample']['collection_location'] = term_to_uri_dict[VALUE_text]
                elif VALUE_text.lower() not in ['na', 'not applicable']:
                    missing_value_list.append('\t'.join([accession, 'geo_loc_name', VALUE_text]))
            #else:
            #    if TAG_text not in ['lat_lon', 'host_disease', 'BioSampleModel', 'passage_history']:
            #        print(accession, TAG_text, VALUE_text)


    taxon_id = SAMPLE.find('SAMPLE_NAME').find('TAXON_ID').text
    info_for_yaml_dict['virus']['virus_species'] = "http://purl.obolibrary.org/obo/NCBITaxon_"+taxon_id


    EXPERIMENT = EXPERIMENT_PACKAGE.find('EXPERIMENT')
    INSTRUMENT_MODEL = [x.text for x in EXPERIMENT.find('PLATFORM').iter('INSTRUMENT_MODEL')][0]

    if INSTRUMENT_MODEL.lower() != 'unspecified':
        if INSTRUMENT_MODEL in term_to_uri_dict:
            info_for_yaml_dict['technology']['sample_sequencing_technology'] = [term_to_uri_dict[INSTRUMENT_MODEL]]
        else:
            missing_value_list.append('\t'.join([accession, 'sample_sequencing_technology', INSTRUMENT_MODEL]))
    #else:
    #    print(accession, 'Missing INSTRUMENT_MODEL', info_for_yaml_dict)
    LIBRARY_DESCRIPTOR = EXPERIMENT.find('DESIGN').find('LIBRARY_DESCRIPTOR')
    if LIBRARY_DESCRIPTOR.text not in ['OTHER']:
        info_for_yaml_dict['technology']['additional_technology_information'] = 'LIBRARY_STRATEGY: {};'.format(LIBRARY_DESCRIPTOR.find('LIBRARY_STRATEGY').text)

    SUBMISSION = EXPERIMENT_PACKAGE.find('SUBMISSION')
    info_for_yaml_dict['submitter']['submitter_sample_id'] = SUBMISSION.attrib['accession']

    if SUBMISSION.attrib['lab_name'].lower() not in ['na']:
        info_for_yaml_dict['submitter']['originating_lab'] = SUBMISSION.attrib['lab_name']

    STUDY = EXPERIMENT_PACKAGE.find('STUDY')
    info_for_yaml_dict['submitter']['publication'] = STUDY.attrib['alias']


    Organization = EXPERIMENT_PACKAGE.find('Organization')
    Organization_Name = Organization.find('Name')
    info_for_yaml_dict['submitter']['authors'] = [Organization_Name.text]

    Organization_Contact = Organization.find('Contact')
    if Organization_Contact is not None:
        Organization_Contact_Name = Organization_Contact.find('Name')
        info_for_yaml_dict['submitter']['submitter_name'] = [Organization_Contact_Name.find('First').text + ' ' + Organization_Contact_Name.find('Last').text]
        info_for_yaml_dict['submitter']['additional_submitter_information'] = Organization_Contact.attrib['email']

        Organization_Concact_Address = Organization_Contact.find('Address')
        if Organization_Concact_Address is not None:
            info_for_yaml_dict['submitter']['submitter_address'] = '; '.join([x.text for x in Organization_Concact_Address] + ['Postal code ' + Organization_Concact_Address.attrib['postal_code']])

    Organization_Address = Organization.find('Address')
    if Organization_Address is not None:
        info_for_yaml_dict['submitter']['lab_address'] = '; '.join([x.text for x in Organization_Address] + ['Postal code ' + Organization_Address.attrib['postal_code']])

    if 'collection_date' not in info_for_yaml_dict['sample']:
        info_for_yaml_dict['sample']['collection_date'] = '1970-01-01'
        info_for_yaml_dict['sample']['additional_collection_information'] = "The real 'collection_date' is missing"


    # Check if mandatory fields are missing
    if 'sample_sequencing_technology' not in info_for_yaml_dict['technology']:
        # print(accession_version, ' - technology not found')
        if accession not in not_created_accession_dict:
            not_created_accession_dict[accession] = []
        not_created_accession_dict[accession].append('sample_sequencing_technology not found')

    if 'collection_location' not in info_for_yaml_dict['sample']:
        if accession not in not_created_accession_dict:
            not_created_accession_dict[accession] = []
        not_created_accession_dict[accession].append('collection_location not found')

    if 'collection_date' not in info_for_yaml_dict['sample']:
        if accession not in not_created_accession_dict:
            not_created_accession_dict[accession] = []
        not_created_accession_dict[accession].append('collection_date not found')

    if 'authors' not in info_for_yaml_dict['submitter']:
        if accession not in not_created_accession_dict:
            not_created_accession_dict[accession] = []
        not_created_accession_dict[accession].append('authors not found')

    if 'host_species' not in info_for_yaml_dict['host']:
        if accession not in not_created_accession_dict:
            not_created_accession_dict[accession] = []
        not_created_accession_dict[accession].append('host_species not found')

    if accession not in not_created_accession_dict:
        num_yaml_created += 1

        with open(os.path.join(dir_yaml, '{}.yaml'.format(accession)), 'w') as fw:
            json.dump(info_for_yaml_dict, fw, indent=2)

if len(missing_value_list) > 0:
    path_missing_terms_tsv = 'missing_terms.sra.tsv'
    print('Written missing terms in {}'.format(path_missing_terms_tsv))
    with open(path_missing_terms_tsv, 'w') as fw:
        fw.write('\n'.join(missing_value_list))

if len(not_created_accession_dict) > 0:
    path_not_created_accession_tsv = 'not_created_accession.sra.tsv'
    print('Written not created accession in {}'.format(path_not_created_accession_tsv))
    with open(path_not_created_accession_tsv, 'w') as fw:
        fw.write('\n'.join(['\t'.join([accession_version, ','.join(missing_info_list)]) for accession_version, missing_info_list in not_created_accession_dict.items()]))

print('Num. YAML files created: {}'.format(num_yaml_created))
