# Genbank XML parser

import os
import sys
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
"""

def get_metadata(id, gb):
    return True,None

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

min_len_to_count = 15000
num_seq_with_len_ge_X_bp = 0

missing_value_list = []
not_created_accession_dict = {}
accession_with_errors_list = []
if None:

    tree = ET.parse(path_metadata_xxx_xml)
    GBSet = tree.getroot()

    for GBSeq in GBSet:
        accession_version = GBSeq.find('GBSeq_accession-version').text

        try:
            info_for_yaml_dict = {
                'id': 'placeholder',
                'host': {},
                'sample': {},
                'virus': {},
                'technology': {},
                'submitter': {}
            }

            sample['sample_id'] = accession_version
            sample['source_database_accession'] = ["http://identifiers.org/insdc/"+accession_version+"#sequence"] #accession is turned into resolvable URL/URI now

            # submitter info
            GBSeq_references = GBSeq.find('GBSeq_references')
            if GBSeq_references is not None:
                author_list = ["{}".format(x.text) for x in GBSeq_references.iter('GBAuthor')]
                if len(author_list) > 0:
                    submitter['authors'] = author_list

                GBReference = GBSeq_references.find('GBReference')
                if GBReference is not None:
                    GBReference_journal = GBReference.find('GBReference_journal')

                    if GBReference_journal is not None and GBReference_journal.text != 'Unpublished':
                        if 'Submitted' in GBReference_journal.text:
                            submitter['submitter_name'] = ["{}".format(GBReference_journal.text.split(') ')[1].split(',')[0].strip())]
                            submitter['submitter_address'] = ','.join(GBReference_journal.text.split(') ')[1].split(',')[1:]).strip()
                        else:
                            submitter['additional_submitter_information'] = GBReference_journal.text

            # This script download and prepare data and metadata for assemblies samples
            technology['assembly_method'] = 'http://purl.obolibrary.org/obo/GENEPIO_0001628'

            GBSeq_comment = GBSeq.find('GBSeq_comment')
            if GBSeq_comment is not None and 'Assembly-Data' in GBSeq_comment.text:
                prefix_split_string = '##Genome-Assembly' if GBSeq_comment.text.startswith('##Genome-') else '##Assembly'

                GBSeq_comment_text = GBSeq_comment.text.split(
                    '{}-Data-START## ; '.format(prefix_split_string)
                )[1].split(' ; {}-Data-END##'.format(prefix_split_string))[0]

                for info_to_check, field_in_yaml in zip(
                    ['Assembly Method', 'Coverage', 'Sequencing Technology'],
                    ['alignment_protocol', 'sequencing_coverage', 'sample_sequencing_technology']
                ):
                    if info_to_check in GBSeq_comment_text:
                        tech_info_to_parse = GBSeq_comment_text.split('{} :: '.format(info_to_check))[1].split(' ;')[0]

                        if field_in_yaml == 'sequencing_coverage':
                            # A regular expression would be better!
                            try:
                                technology[field_in_yaml] = [
                                    float(tech_info_to_parse.replace('(average)', '').replace("reads/nt", '').
                                          replace('(average for 6 sequences)', '').replace(',', '.').strip(' xX>'))
                                ]
                            except ValueError:
                                print(accession_version, "Couldn't make sense of Coverage '%s'" % tech_info_to_parse)
                                pass
                        elif field_in_yaml == 'sample_sequencing_technology':
                            new_seq_tec_list = []
                            for seq_tec in tech_info_to_parse.split(';'):
                                seq_tec = seq_tec.strip()
                                if seq_tec in field_to_term_to_uri_dict['ncbi_sequencing_technology']:
                                    seq_tec = field_to_term_to_uri_dict['ncbi_sequencing_technology'][seq_tec]
                                    new_seq_tec_list.append(seq_tec)
                                else:
                                    missing_value_list.append('\t'.join([accession_version, 'sample_sequencing_technology', seq_tec]))

                            if len(new_seq_tec_list) > 0:
                                technology['sample_sequencing_technology'] = [x for x in new_seq_tec_list]
                        else:
                            technology[field_in_yaml] = tech_info_to_parse


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
                        if GBQualifier_value_text in field_to_term_to_uri_dict['ncbi_host_species']:
                            # Cases like 'Felis catus; Domestic Shorthair'
                            host['host_species'] = field_to_term_to_uri_dict['ncbi_host_species'][GBQualifier_value_text]
                        else:
                            GBQualifier_value_text_list = GBQualifier_value_text.split('; ')

                            if GBQualifier_value_text_list[0] in field_to_term_to_uri_dict['ncbi_host_species']:
                                host['host_species'] = field_to_term_to_uri_dict['ncbi_host_species'][GBQualifier_value_text_list[0]]
                            elif GBQualifier_value_text_list[0] and ('MT215193' in accession_version or 'MT270814' in accession_version):
                                # Information checked manually from NCBI Virus
                                host['host_species'] = field_to_term_to_uri_dict['ncbi_host_species']['Canis lupus familiaris']
                            else:
                                missing_value_list.append('\t'.join([accession_version, 'host_species', GBQualifier_value_text_list[0]]))

                            # Possible cases:
                            # - Homo sapiens                        --> ['Homo sapiens']
                            # - Homo sapiens; female                --> ['Homo sapiens', 'female']
                            # - Homo sapiens; female 63             --> ['Homo sapiens', 'female 63']
                            # - Homo sapiens; female; age 40        --> ['Homo sapiens', 'female', 'age 40']
                            # - Homo sapiens; gender: F; age: 61    --> ['Homo sapiens', 'gender: F', 'age: 61']
                            # - Homo sapiens; gender: M; age: 68    --> ['Homo sapiens', 'gender: M', 'age: 68']
                            # - Homo sapiens; hospitalized patient  --> ['Homo sapiens', 'hospitalized patient']
                            # - Homo sapiens; male                  --> ['Homo sapiens', 'male']
                            # - Homo sapiens; male; 63              --> ['Homo sapiens', 'male', '63']
                            # - Homo sapiens; male; age 29          --> ['Homo sapiens', 'male', 'age 29']
                            # - Homo sapiens; symptomatic           --> ['Homo sapiens', 'symptomatic']
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
                                    host['host_sex'] = "http://purl.obolibrary.org/obo/PATO_0000384" if host_sex == 'male' else "http://purl.obolibrary.org/obo/PATO_0000383"
                                elif GBQualifier_value_text_list[1] in field_to_term_to_uri_dict['ncbi_host_health_status']:
                                    host['host_health_status'] = field_to_term_to_uri_dict['ncbi_host_health_status'][GBQualifier_value_text_list[1]]
                                else:
                                    missing_value_list.append('\t'.join([accession_version, 'host_sex or host_health_status', GBQualifier_value_text_list[1]]))

                                # Host age
                                host_age = -1
                                if len(GBQualifier_value_text_list[1].split(' ')) > 1 and is_integer(GBQualifier_value_text_list[1].split(' ')[-1]):
                                    host_age = int(GBQualifier_value_text_list[1].split(' ')[-1])
                                elif len(GBQualifier_value_text_list) > 2 and is_integer(GBQualifier_value_text_list[2].split(' ')[-1]):
                                    host_age = int(GBQualifier_value_text_list[2].split(' ')[-1])

                                if host_age >= 0 and host_age < 110:
                                    host['host_age'] = host_age
                                    host['host_age_unit'] = 'http://purl.obolibrary.org/obo/UO_0000036'
                                elif len(GBQualifier_value_text_list) > 2:
                                    missing_value_list.append('\t'.join([accession_version, 'host_age', GBQualifier_value_text_list[2]]))
                    elif GBQualifier_name_text == 'collected_by':
                        if any([x in GBQualifier_value_text.lower() for x in ['institute', 'hospital', 'city', 'center']]):
                            sample['collecting_institution'] = GBQualifier_value_text
                        else:
                            sample['collector_name'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'isolation_source':
                        if GBQualifier_value_text.upper() in field_to_term_to_uri_dict['ncbi_speciesman_source']:
                            GBQualifier_value_text = GBQualifier_value_text.upper()  # For example, in case of 'usa: wa'

                        # Little cleaning
                        GBQualifier_value_text = GBQualifier_value_text.strip("/'")

                        if GBQualifier_value_text in field_to_term_to_uri_dict['ncbi_speciesman_source']:
                            sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source'][GBQualifier_value_text]]
                        else:
                            if GBQualifier_value_text.lower() in ['np/op', 'np-op', 'np/op swab', 'np/np swab', 'nasopharyngeal and oropharyngeal swab', 'nasopharyngeal/oropharyngeal swab', 'combined nasopharyngeal and oropharyngeal swab', 'naso and/or oropharyngeal swab']:
                                sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasopharyngeal swab'], field_to_term_to_uri_dict['ncbi_speciesman_source']['oropharyngeal swab']]
                            elif GBQualifier_value_text.lower() in ['nasopharyngeal swab/throat swab', 'nasopharyngeal/throat swab', 'nasopharyngeal swab and throat swab', 'nasal swab and throat swab', 'nasopharyngeal aspirate/throat swab', 'Nasopharyngeal/Throat']:
                                sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasopharyngeal swab'], field_to_term_to_uri_dict['ncbi_speciesman_source']['throat swab']]
                            elif GBQualifier_value_text.lower() in ['nasopharyngeal aspirate & throat swab', 'nasopharyngeal aspirate and throat swab']:
                                sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasopharyngeal aspirate'], field_to_term_to_uri_dict['ncbi_speciesman_source']['throat swab']]
                            elif GBQualifier_value_text.lower() in ['nasal swab and throat swab']:
                                sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasal swab'], field_to_term_to_uri_dict['ncbi_speciesman_source']['throat swab']]
                            elif GBQualifier_value_text.lower() in ['nasal-swab and oro-pharyngeal swab']:
                                sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasal swab'], field_to_term_to_uri_dict['ncbi_speciesman_source']['oropharyngeal swab']]
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

                            if 'additional_collection_information' in sample:
                                sample['additional_collection_information'] += "; The 'collection_date' is estimated (the original date was: {})".format(GBQualifier_value_text)
                            else:
                                sample['additional_collection_information'] = "The 'collection_date' is estimated (the original date was: {})".format(GBQualifier_value_text)
                        elif len(GBQualifier_value_text.split('-')) == 2:
                            date_to_write = parse(GBQualifier_value_text).strftime('%Y-%m') + '-15'

                            if 'additional_collection_information' in sample:
                                sample['additional_collection_information'] += "; The 'collection_date' is estimated (the original date was: {})".format(GBQualifier_value_text)
                            else:
                                sample['additional_collection_information'] = "The 'collection_date' is estimated (the original date was: {})".format(GBQualifier_value_text)
                        elif len(GBQualifier_value_text.split('-')) == 3:
                            GBQualifier_value_text_list = GBQualifier_value_text.split('-')

                            if GBQualifier_value_text_list[1].isalpha():
                                date_to_write = parse(GBQualifier_value_text).strftime('%Y-%m-%d')

                        sample['collection_date'] = date_to_write
                    elif GBQualifier_name_text in ['lat_lon', 'country']:
                        if GBQualifier_name_text == 'country' and ': ' in GBQualifier_value_text:
                            GBQualifier_value_text = GBQualifier_value_text.replace(': ', ':')

                        if GBQualifier_value_text in field_to_term_to_uri_dict['ncbi_countries']:
                            sample['collection_location'] = field_to_term_to_uri_dict['ncbi_countries'][GBQualifier_value_text]
                        else:
                            missing_value_list.append('\t'.join([accession_version, GBQualifier_name_text, GBQualifier_value_text]))
                    elif GBQualifier_name_text == 'note':
                        if 'additional_collection_information' in sample:
                            sample['additional_collection_information'] += '; ' + GBQualifier_value_text
                        else:
                            sample['additional_collection_information'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'isolate':
                        virus['virus_strain'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'db_xref':
                        virus['virus_species'] = "http://purl.obolibrary.org/obo/NCBITaxon_"+GBQualifier_value_text.split('taxon:')[1]

            # Check if mandatory fields are missing
            if 'sample_sequencing_technology' not in technology:
                # print(accession_version, ' - technology not found')
                if accession_version not in not_created_accession_dict:
                    not_created_accession_dict[accession_version] = []
                not_created_accession_dict[accession_version].append('sample_sequencing_technology not found')

            if 'collection_location' not in sample:
                if accession_version not in not_created_accession_dict:
                    not_created_accession_dict[accession_version] = []
                not_created_accession_dict[accession_version].append('collection_location not found')

            if 'collection_date' not in sample:
                if accession_version not in not_created_accession_dict:
                    not_created_accession_dict[accession_version] = []
                not_created_accession_dict[accession_version].append('collection_date not found')
            else:
                year, month, day = [int(x) for x in sample['collection_date'].split('-')]

                collection_date_in_yaml = datetime(year, month, day)
                if collection_date_in_yaml < min_acceptable_collection_date:
                    if accession_version not in not_created_accession_dict:
                        not_created_accession_dict[accession_version] = []
                    not_created_accession_dict[accession_version].append('collection_date too early')

            if 'authors' not in submitter:
                if accession_version not in not_created_accession_dict:
                    not_created_accession_dict[accession_version] = []
                not_created_accession_dict[accession_version].append('authors not found')

            if 'host_species' not in host:
                if accession_version not in not_created_accession_dict:
                    not_created_accession_dict[accession_version] = []
                not_created_accession_dict[accession_version].append('host_species not found')

            if len(GBSeq_sequence.text) < min_len_to_count:
                if accession_version not in not_created_accession_dict:
                    not_created_accession_dict[accession_version] = []
                not_created_accession_dict[accession_version].append('sequence shorter than {} bp'.format(min_len_to_count))

            if accession_version not in not_created_accession_dict:
                num_seq_with_len_ge_X_bp += 1

                # with open(os.path.join(dir_fasta_and_yaml, '{}.fasta'.format(accession_version)), 'w') as fw:
                #    fw.write('>{}\n{}'.format(accession_version, GBSeq_sequence.text.upper()))

                with open(os.path.join(dir_fasta_and_yaml, '{}.yaml'.format(accession_version)), 'w') as fw:
                    json.dump(info_for_yaml_dict, fw, indent=2)
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

if len(not_created_accession_dict) > 0:
    path_not_created_accession_tsv = 'not_created_accession.genbank.tsv'
    print('Written not created accession in {}'.format(path_not_created_accession_tsv))
    with open(path_not_created_accession_tsv, 'w') as fw:
        fw.write('\n'.join(['\t'.join([accession_version, ','.join(missing_info_list)]) for accession_version, missing_info_list in not_created_accession_dict.items()]))

print('Num. new sequences with length >= {} bp: {}'.format(min_len_to_count, num_seq_with_len_ge_X_bp))
