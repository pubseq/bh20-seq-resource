#!/usr/bin/env python3

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--skip-request', action='store_true', help='skip metadata and sequence request', required=False)
parser.add_argument('--only-missing-id', action='store_true', help='download only missing id', required=False)
parser.add_argument('--dict-ontology', type=str, help='where is the ontology',
                    default='../dict_ontology_standardization/',required=False)
args = parser.parse_args()

from Bio import Entrez
Entrez.email = 'another_email@gmail.com'

import xml.etree.ElementTree as ET
import json
import os
import requests
import sys

from datetime import date
from dateutil.parser import parse

num_ids_for_request = 100

dir_metadata = 'metadata_from_nuccore'
dir_fasta_and_yaml = 'fasta_and_yaml'
dir_dict_ontology_standardization = args.dict_ontology

today_date = date.today().strftime("%Y.%m.%d")
path_ncbi_virus_accession = 'sequences.{}.acc'.format(today_date)

def is_integer(string_to_check):
    try:
        int(string_to_check)
        return True
    except ValueError:
        return False

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

if os.path.exists(dir_metadata):
    print("The directory '{}' already exists.".format(dir_metadata))

    if not args.skip_request:
        print("\tTo start the request, delete the directory '{}' or specify --skip-request.".format(dir_metadata))
        sys.exit(-1)


accession_already_downloaded_set = []

if os.path.exists(dir_fasta_and_yaml):
    print("The directory '{}' already exists.".format(dir_fasta_and_yaml))
    if not args.only_missing_id:
        print("To start the download, delete the directory '{}' or specify --only-missing-id.".format(dir_fasta_and_yaml))
        sys.exit(-1)

    accession_already_downloaded_set = set([x.split('.yaml')[0].split('.')[0] for x in os.listdir(dir_fasta_and_yaml) if x.endswith('.yaml')])
    print('There are {} accession already downloaded.'.format(len(accession_already_downloaded_set)))


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
        tmp_list = [x for x in tmp_list if x[:2] not in ['NM', 'NR', 'NP', 'XM', 'XR', 'XP', 'WP']]

        # Remove the version in the id
        tmp_list = [x.split('.')[0] for x in tmp_list]

        #tmp_list = tmp_list[0:2] # restricting to small run
        new_ids_set = set([x.split('.')[0] for x in tmp_list])
        new_ids = len(new_ids_set.difference(id_set))
        id_set.update(new_ids_set)

        print('Term:', term, '-->', new_ids, 'new IDs from', len(tmp_list), '---> Total unique IDs:', len(id_set))

    if not os.path.exists(path_ncbi_virus_accession):
        r = requests.get('https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?q=*:*&fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=VirusLineageId_ss:(2697049)&cmd=download&sort=SourceDB_s%20desc,CreateDate_dt%20desc,id%20asc&dlfmt=acc&fl=id')
        with open(path_ncbi_virus_accession, 'w') as fw:
            fw.write(r.text)

    with open(path_ncbi_virus_accession) as f:
        tmp_list = [line.strip('\n') for line in f]

    new_ids = len(set(tmp_list).difference(id_set))
    id_set.update(tmp_list)

    print('DB: NCBI Virus', today_date, '-->', new_ids, 'new IDs from', len(tmp_list), '---> Total unique IDs:', len(id_set))

    if len(accession_already_downloaded_set) > 0:
        id_set = id_set.difference(accession_already_downloaded_set)
        print('There are {} missing IDs to download.'.format(len(id_set)))

    os.makedirs(dir_metadata)
    for i, id_x_list in enumerate(chunks(list(id_set), num_ids_for_request)):
        path_metadata_xxx_xml = os.path.join(dir_metadata, 'metadata_{}.xml'.format(i))
        print('Requesting {} ids --> {}'.format(len(id_x_list), path_metadata_xxx_xml))

        with open(path_metadata_xxx_xml, 'w') as fw:
            fw.write(
                Entrez.efetch(db='nuccore', id=id_x_list, retmode='xml').read()
            )


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

if not os.path.exists(dir_fasta_and_yaml):
    os.makedirs(dir_fasta_and_yaml)

min_len_to_count = 27500
num_seq_with_len_ge_X_bp = 0

missing_value_list = []
not_created_accession_list = []
accession_with_errors_list = []

for path_metadata_xxx_xml in [os.path.join(dir_metadata, name_metadata_xxx_xml) for name_metadata_xxx_xml in os.listdir(dir_metadata) if name_metadata_xxx_xml.endswith('.xml')]:
    tree = ET.parse(path_metadata_xxx_xml)
    GBSet = tree.getroot()

    for GBSeq in GBSet:
        accession_version = GBSeq.find('GBSeq_accession-version').text

        GBSeq_sequence = GBSeq.find('GBSeq_sequence')
        if GBSeq_sequence is None:
            print(accession_version, ' - sequence not found')
            continue

        try:
            #print(path_metadata_xxx_xml, accession_version)

            # A general default-empty yaml could be read from the definitive one
            info_for_yaml_dict = {
                'id': 'placeholder',
                'host': {},
                'sample': {},
                'virus': {},
                'technology': {},
                'submitter': {}
            }


            info_for_yaml_dict['sample']['sample_id'] = accession_version
            info_for_yaml_dict['sample']['source_database_accession'] = ["http://identifiers.org/insdc/"+accession_version+"#sequence"] #accession is turned into resolvable URL/URI now


            # submitter info
            GBSeq_references = GBSeq.find('GBSeq_references')
            if GBSeq_references is not None:
                author_list = ["{}".format(x.text) for x in GBSeq_references.iter('GBAuthor')]
                if len(author_list) > 0:
                    info_for_yaml_dict['submitter']['authors'] = author_list

                GBReference = GBSeq_references.find('GBReference')
                if GBReference is not None:
                    GBReference_journal = GBReference.find('GBReference_journal')

                    if GBReference_journal is not None and GBReference_journal.text != 'Unpublished':
                        if 'Submitted' in GBReference_journal.text:
                            info_for_yaml_dict['submitter']['submitter_name'] = ["{}".format(GBReference_journal.text.split(') ')[1].split(',')[0].strip())]
                            info_for_yaml_dict['submitter']['submitter_address'] = ','.join(GBReference_journal.text.split(') ')[1].split(',')[1:]).strip()
                        else:
                            info_for_yaml_dict['submitter']['additional_submitter_information'] = GBReference_journal.text


            GBSeq_comment = GBSeq.find('GBSeq_comment')
            if GBSeq_comment is not None and 'Assembly-Data' in GBSeq_comment.text:
                prefix_split_string = '##Genome-Assembly' if GBSeq_comment.text.startswith('##Genome-') else '##Assembly'

                GBSeq_comment_text = GBSeq_comment.text.split(
                    '{}-Data-START## ; '.format(prefix_split_string)
                )[1].split(' ; {}-Data-END##'.format(prefix_split_string))[0]

                for info_to_check, field_in_yaml in zip(
                    ['Assembly Method', 'Coverage', 'Sequencing Technology'],
                    ['sequence_assembly_method', 'sequencing_coverage', 'sample_sequencing_technology']
                ):
                    if info_to_check in GBSeq_comment_text:
                        tech_info_to_parse = GBSeq_comment_text.split('{} :: '.format(info_to_check))[1].split(' ;')[0]

                        if field_in_yaml == 'sequencing_coverage':
                            # A regular expression would be better!
                            try:
                                info_for_yaml_dict['technology'][field_in_yaml] = [
                                    float(tech_info_to_parse.strip('(average)').strip("reads/nt").strip('(average for 6 sequences)').replace(',', '.').strip(' xX>'))
                                ]
                            except ValueError:
                                print(accession_version, "Couldn't make sense of Coverage '%s'" % tech_info_to_parse)
                                pass
                        elif field_in_yaml == 'sample_sequencing_technology':
                            new_seq_tec_list = []
                            for seq_tec in tech_info_to_parse.split(';'):
                                seq_tec = seq_tec.strip()
                                if seq_tec in term_to_uri_dict:
                                    seq_tec = term_to_uri_dict[seq_tec]
                                    new_seq_tec_list.append(seq_tec)
                                else:
                                    missing_value_list.append('\t'.join([accession_version, 'sample_sequencing_technology', seq_tec]))

                            if len(new_seq_tec_list) > 0:
                                info_for_yaml_dict['technology']['sample_sequencing_technology'] = [x for x in new_seq_tec_list]
                        else:
                            info_for_yaml_dict['technology'][field_in_yaml] = tech_info_to_parse


            for GBFeature in GBSeq.iter('GBFeature'):
                if GBFeature.find('GBFeature_key').text != 'source':
                    continue

                for GBQualifier in GBFeature.iter('GBQualifier'):
                    GBQualifier_value = GBQualifier.find('GBQualifier_value')
                    if GBQualifier_value is None:
                        continue
                    GBQualifier_value_text = GBQualifier_value.text

                    GBQualifier_name_text = GBQualifier.find('GBQualifier_name').text

                    if GBQualifier_name_text == 'host':
                        GBQualifier_value_text = GBQualifier_value_text.split(';')[0] # For case like Homo sapiens;sex:female
                        if GBQualifier_value_text in term_to_uri_dict:
                            # Cases like 'Felis catus; Domestic Shorthair'
                            info_for_yaml_dict['host']['host_species'] = term_to_uri_dict[GBQualifier_value_text]
                        else:
                            GBQualifier_value_text_list = GBQualifier_value_text.split('; ')

                            if GBQualifier_value_text_list[0] in term_to_uri_dict:
                                info_for_yaml_dict['host']['host_species'] = term_to_uri_dict[GBQualifier_value_text_list[0]]
                            elif GBQualifier_value_text_list[0] and ('MT215193' in accession_version or 'MT270814' in accession_version):
                                # Information checked manually from NCBI Virus
                                info_for_yaml_dict['host']['host_species'] = term_to_uri_dict['Canis lupus familiaris']
                            else:
                                missing_value_list.append('\t'.join([accession_version, 'host_species', GBQualifier_value_text_list[0]]))

                            # Possible cases:
                            # - Homo sapiens						--> ['Homo sapiens']
                            # - Homo sapiens; female				--> ['Homo sapiens', 'female']
                            # - Homo sapiens; female 63				--> ['Homo sapiens', 'female 63']
                            # - Homo sapiens; female; age 40		--> ['Homo sapiens', 'female', 'age 40']
                            # - Homo sapiens; gender: F; age: 61	--> ['Homo sapiens', 'gender: F', 'age: 61']
                            # - Homo sapiens; gender: M; age: 68	--> ['Homo sapiens', 'gender: M', 'age: 68']
                            # - Homo sapiens; hospitalized patient	--> ['Homo sapiens', 'hospitalized patient']
                            # - Homo sapiens; male					--> ['Homo sapiens', 'male']
                            # - Homo sapiens; male; 63				--> ['Homo sapiens', 'male', '63']
                            # - Homo sapiens; male; age 29			--> ['Homo sapiens', 'male', 'age 29']
                            # - Homo sapiens; symptomatic			--> ['Homo sapiens', 'symptomatic']
                            if len(GBQualifier_value_text_list) > 1:
                                host_sex = ''
                                if 'female' in GBQualifier_value_text_list[1]:
                                    host_sex = 'female'
                                elif 'male' in GBQualifier_value_text_list[1]:
                                    host_sex = 'male'
                                elif 'gender' in GBQualifier_value_text_list[1]:
                                    host_sex_one_lecter = GBQualifier_value_text_list[1].split(':')[-1].strip()
                                    if host_sex_one_lecter in ['F', 'M']:
                                        host_sex = 'female' if host_sex_one_lecter == 'F' else 'male'

                                if host_sex in ['male', 'female']:
                                    info_for_yaml_dict['host']['host_sex'] = "http://purl.obolibrary.org/obo/PATO_0000384" if host_sex == 'male' else "http://purl.obolibrary.org/obo/PATO_0000383"
                                elif GBQualifier_value_text_list[1] in term_to_uri_dict:
                                    info_for_yaml_dict['host']['host_health_status'] = term_to_uri_dict[GBQualifier_value_text_list[1]]
                                else:
                                    missing_value_list.append('\t'.join([accession_version, 'host_sex or host_health_status', GBQualifier_value_text_list[1]]))

                                # Host age
                                host_age = -1
                                if len(GBQualifier_value_text_list[1].split(' ')) > 1 and is_integer(GBQualifier_value_text_list[1].split(' ')[-1]):
                                    host_age = int(GBQualifier_value_text_list[1].split(' ')[-1])
                                elif len(GBQualifier_value_text_list) > 2 and is_integer(GBQualifier_value_text_list[2].split(' ')[-1]):
                                    host_age = int(GBQualifier_value_text_list[2].split(' ')[-1])

                                if host_age > -1:
                                    info_for_yaml_dict['host']['host_age'] = host_age
                                    info_for_yaml_dict['host']['host_age_unit'] = 'http://purl.obolibrary.org/obo/UO_0000036'
                                elif len(GBQualifier_value_text_list) > 2:
                                    missing_value_list.append('\t'.join([accession_version, 'host_age', GBQualifier_value_text_list[2]]))
                    elif GBQualifier_name_text == 'collected_by':
                        if any([x in GBQualifier_value_text.lower() for x in ['institute', 'hospital', 'city', 'center']]):
                            info_for_yaml_dict['sample']['collecting_institution'] = GBQualifier_value_text
                        else:
                            info_for_yaml_dict['sample']['collector_name'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'isolation_source':
                        if GBQualifier_value_text.upper() in term_to_uri_dict:
                            GBQualifier_value_text = GBQualifier_value_text.upper() # For example, in case of 'usa: wa'

                        # Little cleaning
                        GBQualifier_value_text = GBQualifier_value_text.strip("/'")

                        if GBQualifier_value_text in term_to_uri_dict:
                            info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict[GBQualifier_value_text]]
                        else:
                            if GBQualifier_value_text.lower() in ['np/op', 'np/op swab', 'np/np swab', 'nasopharyngeal and oropharyngeal swab', 'nasopharyngeal/oropharyngeal swab', 'combined nasopharyngeal and oropharyngeal swab', 'naso and/or oropharyngeal swab']:
                                info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasopharyngeal swab'], term_to_uri_dict['oropharyngeal swab']]
                            elif GBQualifier_value_text.lower() in ['nasopharyngeal swab/throat swab', 'nasopharyngeal/throat swab', 'nasopharyngeal swab and throat swab', 'nasal swab and throat swab', 'nasopharyngeal aspirate/throat swab', 'Nasopharyngeal/Throat']:
                                info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasopharyngeal swab'], term_to_uri_dict['throat swab']]
                            elif GBQualifier_value_text.lower() in ['nasopharyngeal aspirate & throat swab', 'nasopharyngeal aspirate and throat swab']:
                                info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasopharyngeal aspirate'], term_to_uri_dict['throat swab']]
                            elif GBQualifier_value_text.lower() in ['nasal swab and throat swab']:
                                info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasal swab'], term_to_uri_dict['throat swab']]
                            elif GBQualifier_value_text.lower() in ['nasal-swab and oro-pharyngeal swab']:
                                info_for_yaml_dict['sample']['specimen_source'] = [term_to_uri_dict['nasal swab'], term_to_uri_dict['oropharyngeal swab']]
                            else:
                                missing_value_list.append('\t'.join([accession_version, 'specimen_source', GBQualifier_value_text]))
                    elif GBQualifier_name_text == 'collection_date':
                        # TO_DO: which format we will use?
                        date_to_write = GBQualifier_value_text

                        if len(GBQualifier_value_text.split('-')) == 1:
                            if int(GBQualifier_value_text) < 2020:
                                date_to_write = "{}-12-15".format(GBQualifier_value_text)
                            else:
                                date_to_write = "{}-01-15".format(GBQualifier_value_text)

                            if 'additional_collection_information' in info_for_yaml_dict['sample']:
                                info_for_yaml_dict['sample']['additional_collection_information'] += "; The 'collection_date' is estimated (the original date was: {})".format(GBQualifier_value_text)
                            else:
                                info_for_yaml_dict['sample']['additional_collection_information'] = "The 'collection_date' is estimated (the original date was: {})".format(GBQualifier_value_text)
                        elif len(GBQualifier_value_text.split('-')) == 2:
                            date_to_write += '-15'

                            if 'additional_collection_information' in info_for_yaml_dict['sample']:
                                info_for_yaml_dict['sample']['additional_collection_information'] += "; The 'collection_date' is estimated (the original date was: {})".format(GBQualifier_value_text)
                            else:
                                info_for_yaml_dict['sample']['additional_collection_information'] = "The 'collection_date' is estimated (the original date was: {})".format(GBQualifier_value_text)
                        elif len(GBQualifier_value_text.split('-')) == 3:
                            GBQualifier_value_text_list = GBQualifier_value_text.split('-')

                            if GBQualifier_value_text_list[1].isalpha():
                                date_to_write = parse(GBQualifier_value_text).strftime('%Y-%m-%d')

                        info_for_yaml_dict['sample']['collection_date'] = date_to_write
                    elif GBQualifier_name_text in ['lat_lon', 'country']:
                        if GBQualifier_name_text == 'country' and ': ' in GBQualifier_value_text:
                            GBQualifier_value_text = GBQualifier_value_text.replace(': ', ':')

                        if GBQualifier_value_text in term_to_uri_dict:
                            info_for_yaml_dict['sample']['collection_location'] = term_to_uri_dict[GBQualifier_value_text]
                        else:
                            missing_value_list.append('\t'.join([accession_version, GBQualifier_name_text, GBQualifier_value_text]))
                    elif GBQualifier_name_text == 'note':
                        if 'additional_collection_information' in info_for_yaml_dict['sample']:
                            info_for_yaml_dict['sample']['additional_collection_information'] += '; ' + GBQualifier_value_text
                        else:
                            info_for_yaml_dict['sample']['additional_collection_information'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'isolate':
                        info_for_yaml_dict['virus']['virus_strain'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'db_xref':
                        info_for_yaml_dict['virus']['virus_species'] = "http://purl.obolibrary.org/obo/NCBITaxon_"+GBQualifier_value_text.split('taxon:')[1]


            if 'sample_sequencing_technology' not in info_for_yaml_dict['technology']:
                #print(accession_version, ' - technology not found')
                not_created_accession_list.append([accession_version, 'technology not found'])
                continue

            with open(os.path.join(dir_fasta_and_yaml, '{}.fasta'.format(accession_version)), 'w') as fw:
                fw.write('>{}\n{}'.format(accession_version, GBSeq_sequence.text.upper()))

            with open(os.path.join(dir_fasta_and_yaml, '{}.yaml'.format(accession_version)), 'w') as fw:
                json.dump(info_for_yaml_dict, fw, indent=2)


            if(len(GBSeq_sequence.text) >= min_len_to_count):
                num_seq_with_len_ge_X_bp += 1
        except:
            print("Unexpected error for the ID {}: {}".format(accession_version, sys.exc_info()[0]))
            accession_with_errors_list.append(accession_version)
            continue

if len(missing_value_list) > 0:
    path_missing_terms_tsv = 'missing_terms.genbank.tsv'
    print('Written missing terms in {}'.format(path_missing_terms_tsv))
    with open(path_missing_terms_tsv, 'w') as fw:
        fw.write('\n'.join(missing_value_list))

if len(accession_with_errors_list) > 0:
    path_accession_with_errors_tsv = 'accession_with_errors.genbank.tsv'
    print('Written the accession with errors in {}'.format(path_accession_with_errors_tsv))
    with open(path_accession_with_errors_tsv, 'w') as fw:
        fw.write('\n'.join(accession_with_errors_list))

if len(not_created_accession_list) > 0:
    path_not_created_accession_tsv = 'not_created_accession.genbank.tsv'
    print('Written not created accession in {}'.format(path_not_created_accession_tsv))
    with open(path_not_created_accession_tsv, 'w') as fw:
        fw.write('\n'.join(['\t'.join(x) for x in not_created_accession_list]))

print('Num. new sequences with length >= {} bp: {}'.format(min_len_to_count, num_seq_with_len_ge_X_bp))
