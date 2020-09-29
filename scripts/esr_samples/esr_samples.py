import os
import pandas as pd
from string import Template
from dateutil.parser import parse

import sys

sys.path.append('../')
from utils import check_and_get_ontology_dictionaries

# Metadata in tabular format
path_metadata_xlsx = 'Pathogen.cl.1.0.xlsx'

path_template_yaml = 'template.yaml'
# Removed from the template (for now)
# license:
#    license_type: "http://creativecommons.org/licenses/by/4.0/"
#    title: "SARS-CoV-2 New Zealand"
#    attribution_name: "ESR"
#    attribution_url: "https://www.esr.cri.nz/"


# Read the dictionaries for the ontology
dir_dict_ontology_standardization = '../dict_ontology_standardization/'
field_to_term_to_uri_dict = check_and_get_ontology_dictionaries(dir_dict_ontology_standardization)

dir_output = 'yaml'
suffix = '.consensus'

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

metadata_df = pd.read_excel(path_metadata_xlsx, skiprows=12)

# Maybe not the best pandas-way to do this
for index, row in metadata_df.iterrows():
    # print(row['*sample_name'])

    geo_loc_name = row['*geo_loc_name'].replace(': ', ':')

    if geo_loc_name not in field_to_term_to_uri_dict['ncbi_countries']:
        if geo_loc_name in [
            'New Zealand:Counties Manukau', 'New Zealand:Capital and Coast', 'New Zealand:Southern',
            'New Zealand:Waikato',
            'New Zealand:Lakes', 'New Zealand:Nelson Marlborough', 'New Zealand:South Canterbury',
            'New Zealand:MidCentral',
            'New Zealand:Tairawhiti', 'New Zealand:Hawkes Bay', 'New Zealand:NA', 'New Zealand:Taranaki'
        ]:
            geo_loc_name = 'New Zealand'
        else:
            print(geo_loc_name)
            break

    country = field_to_term_to_uri_dict['ncbi_countries'][geo_loc_name]

    d = {
        'host_species': field_to_term_to_uri_dict['ncbi_host_species'][row['*host']],
        'sample_id': row['*sample_name'],
        'collection_date': parse(row['*collection_date']).strftime('%Y-%m-%d'),
        'collection_location': country,
        'specimen_source': field_to_term_to_uri_dict['ncbi_speciesman_source'][row['*isolation_source']],
        'virus_species': 'http://purl.obolibrary.org/obo/NCBITaxon_2697049',

        'submitter_sample_id': row['bioproject_accession'],
    }

    with open(path_template_yaml) as f:
        src = Template(f.read())

        with open(os.path.join(dir_output, '{}{}.yaml'.format(row['*sample_name'], suffix)), 'w') as fw:
            fw.write(src.substitute(d))

print('{} YAML files created.'.format(len([x for x in os.listdir(dir_output) if x.endswith('.yaml')])))
