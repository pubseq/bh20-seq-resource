from Bio import Entrez
Entrez.email = 'your_email_to_be_polite'

import xml.etree.ElementTree as ET
import yaml
import os

path_ncbi_virus_accession = 'sequences.acc'

date = '20200415'
path_seq_fasta = 'seq_from_nuccore.{}.fasta'.format(date)
path_metadata_xml = 'metadata_from_nuccore.{}.xml'.format(date)

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

    id_set.update([x.split('.')[0] for x in tmp_list])

print(term_list, len(id_set))

with open(path_ncbi_virus_accession) as f:
    tmp_list = [line.strip('\n') for line in f]

print('NCBI Virus', len(tmp_list))
id_set.update(tmp_list)

print(term_list + ['NCBI Virus'], len(id_set))

if not os.path.exists(path_metadata_xml):
    # TO_DO: to check if I already have the records?
    
    with open(path_metadata_xml, 'w') as fw:
        fw.write(
            Entrez.efetch(db='nuccore', id=list(id_set), retmode='xml').read()
        )
        
        
tree = ET.parse(path_metadata_xml)
GBSet = tree.getroot()

species_to_taxid_dict = {
    'Homo sapiens': 9606
}

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
    info_for_yaml_dict['submitter']['authors'] = ';'.join([x.text for x in GBSeq.iter('GBAuthor')])

    
    GBSeq_comment = GBSeq.find('GBSeq_comment')
    if GBSeq_comment is not None and 'Assembly-Data' in GBSeq_comment.text:
        GBSeq_comment_text = GBSeq_comment.text.split('##Assembly-Data-START## ; ')[1].split(' ; ##Assembly-Data-END##')[0]

        for info_to_check, field_in_yaml in zip(
            ['Assembly Method', 'Coverage', 'Sequencing Technology'],
            ['sequence_assembly_method', 'sequencing_coverage', 'sample_sequencing_technology']
        ):
            if info_to_check in GBSeq_comment_text:
                info_for_yaml_dict['technology'][field_in_yaml] = GBSeq_comment_text.split('{} :: '.format(info_to_check))[1].split(' ;')[0]
    
    
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

                info_for_yaml_dict['host']['host_common_name'] = GBQualifier_value_text_list[0]

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
                info_for_yaml_dict['sample']['specimen_source'] = GBQualifier_value_text
            elif GBQualifier_name_text == 'collection_date':
                # TO_DO: which format we will use?
                info_for_yaml_dict['sample']['collection_date'] = GBQualifier_value_text
            elif GBQualifier_name_text in ['lat_lon', 'country']:
                info_for_yaml_dict['sample']['collection_location'] = GBQualifier_value_text
            elif GBQualifier_name_text == 'note':
                info_for_yaml_dict['sample']['additional_collection_information'] = GBQualifier_value_text
            elif GBQualifier_name_text == 'isolate':
                info_for_yaml_dict['virus']['virus_strain'] = GBQualifier_value_text
            elif GBQualifier_name_text == 'db_xref':
                info_for_yaml_dict['virus']['virus_species'] = int(GBQualifier_value_text.split('taxon:')[1])
    
    with open('{}.fasta'.format(accession_version), 'w') as fw:
        fw.write('>{}\n{}'.format(accession_version, GBSeq_sequence.text.upper()))

    with open('{}.yaml'.format(accession_version), 'w') as fw:
        yaml.dump(info_for_yaml_dict, fw, default_flow_style=False)
