# pipeline

```sh
# --- get list of IDs already in PubSeq
./sparql-fetch-ids > pubseq_ids.txt
# --- get list of missing genbank IDs
./genbank-fetch-ids.py --skip pubseq_ids.txt > genbank_ids.txt
# --- fetch XML
python3 update-from-genbank.py --ids genbank_ids.txt --out ~/tmp/genbank
# --- Transform to YAML and FASTA
python3 transform-genbank-xml2yamlfa --out ~/tmp/pubseq file(s)
```

# TODO

- [ ] Add id for GenBank accession - i.e. how can we tell a record is from GenBank
