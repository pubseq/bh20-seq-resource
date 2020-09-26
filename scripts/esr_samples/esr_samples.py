import os
import pandas as pd
from string import Template
from dateutil.parser import parse

path_metadata_xlsx = 'Pathogen.cl.1.0.xlsx'

path_template_yaml = 'template.yaml'
# Removed from the template (for now)
# license:
#    license_type: "http://creativecommons.org/licenses/by/4.0/"
#    title: "SARS-CoV-2 New Zealand"
#    attribution_name: "ESR"
#    attribution_url: "https://www.esr.cri.nz/"

dir_dict_ontology_standardization = '../dict_ontology_standardization/'

dir_output = 'yaml'
suffix = '.consensus'

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

term_to_uri_dict = {}

for path_dict_xxx_csv in [os.path.join(dir_dict_ontology_standardization, name_xxx_csv) for name_xxx_csv in
                          os.listdir(dir_dict_ontology_standardization) if name_xxx_csv.endswith('.csv')]:
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

metadata_df = pd.read_excel(path_metadata_xlsx, skiprows=12)

# Maybe not the best pandas-way to do this
for index, row in metadata_df.iterrows():
    # print(row['*sample_name'])

    geo_loc_name = row['*geo_loc_name'].replace(': ', ':')
    country = ''
    if not geo_loc_name in term_to_uri_dict:
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

    country = term_to_uri_dict[geo_loc_name]

    d = {
        'host_species': term_to_uri_dict[row['*host']],
        'sample_id': row['*sample_name'],
        'collection_date': parse(row['*collection_date']).strftime('%Y-%m-%d'),
        'collection_location': country,
        'specimen_source': term_to_uri_dict[row['*isolation_source']],
        'virus_species': 'http://purl.obolibrary.org/obo/NCBITaxon_2697049',

        'submitter_sample_id': row['bioproject_accession'],
    }

    with open(path_template_yaml) as f:
        src = Template(f.read())

        with open(os.path.join(dir_output, '{}{}.yaml'.format(row['*sample_name'], suffix)), 'w') as fw:
            fw.write(src.substitute(d))

print('{} YAML files created.'.format(len([x for x in os.listdir(dir_output) if x.endswith('.yaml')])))
