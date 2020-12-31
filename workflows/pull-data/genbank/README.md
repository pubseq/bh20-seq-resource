Pipeline:

```sh
# --- get list of IDs already in PubSeq
sparql-fetch-ids > pubseq_ids.txt
# --- fetch XML
update-from-genbank --skip pubseq_ids.txt --max 100 --outdir ~/tmp/genbank
# --- get new IDs
genbank-fetch-ids --dir ~/tmp/pubseq > genbank_ids.txt
# --- loop through IDs (pseudo code)
for id in genbank_ids.txt:
  transform-genbank-xml2yamlfa --dir ~/tmp/genbank id --outdir ~/tmp/pubseq
```
