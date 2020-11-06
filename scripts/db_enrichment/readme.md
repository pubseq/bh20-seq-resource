We have two files in the folder *semantic_enrichment* that are used to enrich the identifier in our triples store with additional information, e.g. human readable labels and semantics (e.g. *What countries are summarizes as a continent*). This describes how to update these two files.

### semantic_enrichment/labels.ttl
Static label about the ontology vocabulary terms we use. This file has to be updated manually. Use the OLS or bioportal to find more information about a used ontology term.

### semantic_enrichment/countries.ttl
File containing information about the countries in our database. Additional information about countries are e.g. the label or GPS coordinates. We enricht the country identifier via wikidata. Please ensure that the .ttl file is valid by e.g. using his online validator (http://ttl.summerofcode.be/).

#### Update process
- What countries (=wikidata identifier) do we have to enrich?
This SPARQL query (http://sparql.genenetwork.org/sparql/) retrieves all countries (ids) from our database that do not have a label yet:


```sparql
SELECT DISTINCT ?geoLocation  WHERE
{
    ?fasta ?x [ <http://purl.obolibrary.org/obo/GAZ_00000448> ?geoLocation] .
    FILTER NOT EXISTS {?geoLocation <http://www.w3.org/2000/01/rdf-schema#label> ?geoLocation_tmp_label}
}
```

[Run query](http://sparql.genenetwork.org/sparql/?default-graph-uri=&query=%0D%0ASELECT+DISTINCT+%3FgeoLocation++WHERE%0D%0A%7B%0D%0A++%3Ffasta+%3Fx+%5B+%3Chttp%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FGAZ_00000448%3E+%3FgeoLocation%5D+.%0D%0A++FILTER+NOT+EXISTS+%7B%3FgeoLocation+%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23label%3E+%3FgeoLocation_tmp_label%7D%0D%0A%7D&format=text%2Fhtml&timeout=0&debug=on&run=+Run+Query+)

- Use the list of identifiers created with the query above as input for the update script *country_enrichment.py*. The script creates a temporary .ttl file in this folder
- Merge the output of the script above manually into the file semantic_enrichment/countries.ttl (TODO: Improve script output so manual intervention no longer needed. Currently there are "double entries" for continents in the output)
