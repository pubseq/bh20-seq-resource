import os
import pandas as pd
from string import Template
from dateutil.parser import parse
import re

import sys

# Metadata in tabular format in a spreadsheet(?!)
xlsx = '../../test/data/10_samples.xlsx'

# Template in a text file
template_yaml = 'template.yaml'

dir_output = 'yaml'

if not os.path.exists(dir_output):
    os.makedirs(dir_output)

table = pd.read_excel(xlsx)

print(table)

for index, row in table.iterrows():
    sample = row['Sample ID']
    print(f"Processing sample {sample}...")

    with open(template_yaml) as f:
      text = Template(f.read())
      with open(os.path.join(dir_output,f"{sample}.yaml"), 'w') as fw:
          sample_id = sample
          sample_name = sample
          collection_date = parse(str(row['Collection Date'])).strftime('%Y-%m-%d')
          locationx = row['City']+", "+row['State']+", USA"
          location = "http://www.wikidata.org/entity/Q16563" # Memphis by default
          map = {
              "Pegram": "http://www.wikidata.org/entity/Q3289517",
              "Alexander": "http://www.wikidata.org/entity/Q79663",
              "Smithville": "http://www.wikidata.org/entity/Q2145339",
              "Nashville": "http://www.wikidata.org/entity/Q23197",
              "Madison": "http://www.wikidata.org/entity/Q494755"
              }

          for name in map:
              p = re.compile(name)
              if p.match(locationx):
                  location = map[name]
                  break

          strain = f"SARS-CoV-2/human/USA/{sample}/2020"
          fw.write(text.substitute(sample_id=sample_id,
                                   sample_name=sample_name,
                                   collection_date=collection_date,
                                   location=location,
                                   locationx=locationx,
                                   strain=strain
                                   ))

          print(f"Run: python3 bh20sequploader/main.py scripts/uthsc_samples/yaml/{sample}.yaml scripts/uthsc_samples/yaml/{sample}.fa")
