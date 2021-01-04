# GenBank

This directory contains the tools to pull and transform
GenBank data.

# Workflows

## Prepare new GenBank data for upload

The following workflow sends GenBank data into PubSeq

```sh
# --- get list of IDs already in PubSeq
../../tools/sparql-fetch-ids > pubseq_ids.txt
# --- get list of missing genbank IDs
python3 genbank-fetch-ids.py --skip pubseq_ids.txt > genbank_ids.txt

# --- fetch XML
python3 update-from-genbank.py --ids genbank_ids.txt --out ~/tmp/genbank

# --- Transform to YAML/JSON and FASTA
python3 transform-genbank-xml2yamlfa.py --out ~/tmp/pubseq file(s)

# --- Normalize data (validation mode)
python3 ../../workflows/tools/normalize-yamlfa.py -s ~/tmp/yamlfa/state.json --species ncbi_host_species.csv --specimen specimen.csv --validate

```

# TODO

- [ ] Add id for GenBank accession - i.e. how can we tell a record is from GenBank
