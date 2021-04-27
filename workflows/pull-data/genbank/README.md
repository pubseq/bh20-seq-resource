# GenBank

This directory contains the tools to pull and transform
GenBank data.

# Workflows

## Prepare new GenBank data for upload to PubSeq

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
python3 transform-genbank-xml2yamlfa.py --out ~/tmp/yamlfa XMLfile(s)

# After the run you can inspect ~/tmp/yamlfa/state.json for issues
```

This script, apart from unpacking XML to JSON, does some minimal
conversions:

1. Checks the dates (and drops on error)
2. Unpacks the host field into species, sex and age (very rare to have
   these listed)
3. Adds and normalizes GenBank database identifiers and plugs in
   indentifiers.org
4. Normalizes virus_species
5. Checks the sequence is longer than 27K base pairs and writes the
   FASTA file

A record may look like

```json
{
    "update_date": "2021-03-16",
    "host": {
        "host_species": "Homo sapiens"
    },
    "sample": {
        "sample_id": "MW751502",
        "database": "https://www.ncbi.nlm.nih.gov/genbank/",
        "source_database_accession": [
            "http://identifiers.org/insdc/MW751502#sequence"
        ],
    },
    "virus": {
        "virus_strain": "SARS-CoV-2/human/USA/MN-MDH-3477/2021",
        "virus_species": "http://purl.obolibrary.org/obo/NCBITaxon_2697049"
    },
    (...)
    "warnings": []
}
```

### First pass normalisation

At this stage we move from GenBank specific to PubSeq generic
normalisation and validation:

```sh
cd ../../pubseq

# --- Normalize data (validation mode) using state.json
python3 normalize-yamlfa.py -s ~/tmp/yamlfa/state.json --species ../../scripts/dict_ontology_standardization/ncbi_host_species.csv --specimen ../../scripts/dict_ontology_standardization/ncbi_speciesman_source.csv --validate
```

This Python script reads the JSON records and normalises the
host_species, specimen_source fields using regular expressions,
so the JSON can look like

```json
{
  (...)
  "host": {
    "host_species": "http://purl.obolibrary.org/obo/NCBITaxon_9606"
  },
  "sample": {
    "sample_id": "MW751502",
    "database": "https://www.ncbi.nlm.nih.gov/genbank/",
    "source_database_accession": [
      "http://identifiers.org/insdc/MW751502#sequence"
    ],
    "collection_location": "USA: Minnesota",
    "collection_date": "2021-02-21",
    "specimen_source": [
      "http://purl.obolibrary.org/obo/NCIT_C13275"
    ]
  },
  (...)
  "warnings": []
}
```

The validation will stop on errors and they are usually worth fixing.
Note the --start switch to skip previously parsed records.

And when all passes, use the writer option with --out to write the
new JSON files to a directory.

### Second pass normalisation (GEO)

The GEO filter normalises geo information as described in
[this blog](../../../doc/blog/covid19-pubseq-location-data.org).

```sh
ruby normalize.rb -s path/state.json
```

### More

Note the latest writeup can be found [here](https://github.com/pubseq/bh20-seq-resource/blob/master/doc/blog/using-covid-19-pubseq-part3.org#example-uploading-bulk-genbank-sequences)
