# GenBank

This directory contains the tools to pull and transform
GenBank data.

# Workflows

## Prepare new GenBank data for upload

The following workflow fetches GenBank data and uploads that into
PubSeq. First (optionally) set up a Guix environment using
this [guix-run](https://github.com/pubseq/bh20-seq-resource/blob/master/workflows/pull-data/genbank/.guix-run). And set the path to point
to the right python3

```sh
export PATH=$GUIX_ENVIRONMENT/bin:$PATH
which python3
/gnu/store/j1c70...-profile/bin/python3
```

With dependencies set this should be a breeze with something like

```sh
# --- get list of IDs already in PubSeq
../../tools/pubseq-fetch-ids > pubseq_ids.txt

# --- get list of missing genbank IDs
python3 genbank-fetch-ids.py --skip pubseq_ids.txt > genbank_ids.txt

# --- fetch XML
python3 update-from-genbank.py --ids genbank_ids.txt --out ~/tmp/genbank

# --- Transform to YAML/JSON and FASTA
python3 transform-genbank-xml2yamlfa.py --out ~/tmp/pubseq file(s)

# --- Normalize data (validation mode)
python3 ../../workflows/tools/normalize-yamlfa.py -s ~/tmp/yamlfa/state.json --species ncbi_host_species.csv --specimen specimen.csv --validate

```

Note the latest writeup can be found [here](https://github.com/pubseq/bh20-seq-resource/blob/master/doc/blog/using-covid-19-pubseq-part3.org#example-uploading-bulk-genbank-sequences)

## Validate GenBank data

To pull the data from PubSeq use the list of pubseq ids generated
above.



# TODO

- [X] Add id for GenBank accession - i.e. how can we tell a record is from GenBank
