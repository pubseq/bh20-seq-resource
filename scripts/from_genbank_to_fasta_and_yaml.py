from Bio import Entrez
Entrez.email = 'another_email@gmail.com'

import xml.etree.ElementTree as ET
import yaml
import os

from datetime import date
today = date.today().strftime("%Y%m%d")

dir_metadata_today = 'metadata_from_nuccore_{}'.format(today)
dir_fasta_and_yaml_today = 'fasta_and_yaml_{}'.format(today)

dir_dict_ontology_standardization = 'dict_ontology_standardization/'

path_ncbi_virus_accession = 'sequences.acc'

# Take all the ids
id_set = set()

term_list = ['SARS-CoV-2', 'SARS-CoV2', 'SARS CoV2', 'SARSCoV2', 'txid2697049[Organism]']
for term in term_list:
    tmp_list = Entrez.read(
        Entrez.esearch(db='nuccore', term=term, idtype='acc', retmax='10000')
    )['IdList']

    # Remove mRNAs, ncRNAs, Proteins, and predicted models (more information here: https://en.wikipedia.org/wiki/RefSeq)
    tmp_list = [x for x in tmp_list if x[:2] not in ['NM', 'NR', 'NP', 'XM', 'XR', 'XP', 'WP']]

    # Remove the version in the id
    tmp_list = [x.split('.')[0] for x in tmp_list]
    
    print(term, len(tmp_list))
    tmp_list=tmp_list
#    tmp_list = tmp_list[0:2] # restricting to small run

    id_set.update([x.split('.')[0] for x in tmp_list])

print(term_list, len(id_set))

with open(path_ncbi_virus_accession) as f:
    tmp_list = [line.strip('\n') for line in f]

print('NCBI Virus', len(tmp_list))
id_set.update(tmp_list)

print(term_list + ['NCBI Virus'], len(id_set))

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
        
num_ids_for_request = 100
if not os.path.exists(dir_metadata_today):
    os.makedirs(dir_metadata_today)
    
    for i, id_x_list in enumerate(chunks(list(id_set), num_ids_for_request)):
        path_metadata_xxx_xml = os.path.join(dir_metadata_today, 'metadata_{}.xml'.format(i))
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
                term = term.strip('"')
            else:    
                term, uri = line.strip('\n').split(',')

            term_to_uri_dict[term] = uri

species_to_taxid_dict = {
    'Homo sapiens': 'http://purl.obolibrary.org/obo/NCBITaxon_9606'
}


if not os.path.exists(dir_fasta_and_yaml_today):
    os.makedirs(dir_fasta_and_yaml_today)

    for path_metadata_xxx_xml in [os.path.join(dir_metadata_today, name_metadata_xxx_xml) for name_metadata_xxx_xml in os.listdir(dir_metadata_today) if name_metadata_xxx_xml.endswith('.xml')]:
        tree = ET.parse(path_metadata_xxx_xml)
        GBSet = tree.getroot()

        for GBSeq in GBSet:
            accession_version = GBSeq.find('GBSeq_accession-version').text

            GBSeq_sequence = GBSeq.find('GBSeq_sequence')
            if GBSeq_sequence is None:
                print(accession_version, ' - sequence not found')
                continue


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
            info_for_yaml_dict['sample']['source_database_accession'] = accession_version
            info_for_yaml_dict['submitter']['authors'] = ';'.join([x.text for x in GBSeq.iter('GBAuthor')])


            GBSeq_comment = GBSeq.find('GBSeq_comment')
            if GBSeq_comment is not None and 'Assembly-Data' in GBSeq_comment.text:
                GBSeq_comment_text = GBSeq_comment.text.split('##Assembly-Data-START## ; ')[1].split(' ; ##Assembly-Data-END##')[0]

                for info_to_check, field_in_yaml in zip(
                    ['Assembly Method', 'Coverage', 'Sequencing Technology'],
                    ['sequence_assembly_method', 'sequencing_coverage', 'sample_sequencing_technology']
                ):
                    if info_to_check in GBSeq_comment_text:
                        tech_info_to_parse = GBSeq_comment_text.split('{} :: '.format(info_to_check))[1].split(' ;')[0]
                        
                        if field_in_yaml == 'sequencing_coverage':
                            # A regular expression would be better!
                            info_for_yaml_dict['technology'][field_in_yaml] = ';'.join(
                                [x.strip('(average)').strip("reads/nt").replace(',', '.').strip(' xX>') for x in tech_info_to_parse.split(';')]
                            )
                        elif field_in_yaml == 'sample_sequencing_technology':
                            new_seq_tec_list = []
                            for seq_tec in tech_info_to_parse.split(';'):
                                seq_tec = seq_tec.strip()
                                if seq_tec in term_to_uri_dict:
                                    seq_tec = term_to_uri_dict[seq_tec]
                                else:
                                    print(accession_version, 'missing technologies:', seq_tec)
 
                                new_seq_tec_list.append(seq_tec)

                            for n, seq_tec in enumerate(new_seq_tec_list):
                                info_for_yaml_dict['technology'][field_in_yaml + ('' if n == 0 else str(n + 1))] = seq_tec
                        else:
                            info_for_yaml_dict['technology'][field_in_yaml] = tech_info_to_parse

                        
                        #term_to_uri_dict

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
                        GBQualifier_value_text_list = GBQualifier_value_text.split('; ')

                        #info_for_yaml_dict['host']['host_common_name'] = GBQualifier_value_text_list[0] # Removed

                        if GBQualifier_value_text_list[0] in species_to_taxid_dict:
                            info_for_yaml_dict['host']['host_species'] = species_to_taxid_dict[GBQualifier_value_text_list[0]]

                        if len(GBQualifier_value_text_list) > 1:
                            if GBQualifier_value_text_list[1] in ['male', 'female']:
                                info_for_yaml_dict['host']['host_sex'] = GBQualifier_value_text_list[1]
                            else:
                                info_for_yaml_dict['host']['host_health_status'] = GBQualifier_value_text_list[1]

                            if 'age' in GBQualifier_value_text:
                                info_for_yaml_dict['host']['host_age'] = int(GBQualifier_value_text_list[2].split('age ')[1])
                                info_for_yaml_dict['host']['host_age_unit'] = 'year'
                    elif GBQualifier_name_text == 'collected_by':
                        if any([x in GBQualifier_value_text.lower() for x in ['institute', 'hospital', 'city', 'center']]):
                            info_for_yaml_dict['sample']['collecting_institution'] = GBQualifier_value_text
                        else:
                            info_for_yaml_dict['sample']['collector_name'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'isolation_source':
                        if GBQualifier_value_text in term_to_uri_dict:
                            info_for_yaml_dict['sample']['specimen_source'] = term_to_uri_dict[GBQualifier_value_text]
                        else:
                            if GBQualifier_value_text in ['NP/OP swab', 'nasopharyngeal and oropharyngeal swab', 'nasopharyngeal/oropharyngeal swab', 'np/np swab', 'np/op']:
                                info_for_yaml_dict['sample']['specimen_source'] = term_to_uri_dict['nasopharyngeal swab']
                                info_for_yaml_dict['sample']['specimen_source2'] = term_to_uri_dict['oropharyngeal swab']
                            else:
                                print(accession_version, 'missing specimen_source:', GBQualifier_value_text)
                    elif GBQualifier_name_text == 'collection_date':
                        # TO_DO: which format we will use?
                        info_for_yaml_dict['sample']['collection_date'] = GBQualifier_value_text
                    elif GBQualifier_name_text in ['lat_lon', 'country']:
                        if GBQualifier_value_text in term_to_uri_dict:
                            GBQualifier_value_text = term_to_uri_dict[GBQualifier_value_text]
                        else:
                            print(accession_version, 'missing {}:'.format(GBQualifier_name_text), GBQualifier_value_text)

                        info_for_yaml_dict['sample']['collection_location'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'note':
                        info_for_yaml_dict['sample']['additional_collection_information'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'isolate':
                        info_for_yaml_dict['virus']['virus_strain'] = GBQualifier_value_text
                    elif GBQualifier_name_text == 'db_xref':
                        info_for_yaml_dict['virus']['virus_species'] = "http://purl.obolibrary.org/obo/NCBITaxon_"+GBQualifier_value_text.split('taxon:')[1]


            #Remove technology key if empty!
            if (info_for_yaml_dict['technology']=={}):
                del info_for_yaml_dict['technology']

            with open(os.path.join(dir_fasta_and_yaml_today, '{}.fasta'.format(accession_version)), 'w') as fw:
                fw.write('>{}\n{}'.format(accession_version, GBSeq_sequence.text.upper()))

            with open(os.path.join(dir_fasta_and_yaml_today, '{}.yaml'.format(accession_version)), 'w') as fw:
                yaml.dump(info_for_yaml_dict, fw, default_flow_style=False)
