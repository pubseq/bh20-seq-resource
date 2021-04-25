#! /bin/bash
#
# This query fetches approx 80K places
#
curl -G https://query.wikidata.org/sparql -H "Accept: text/csv; charset=utf-8" --data-urlencode query="
SELECT DISTINCT ?placename ?place ?country ?coor ?population where {
  ?place wdt:P17 ?country ;
         wdt:P625 ?coor ;
         wdt:P1082 ?population .
  FILTER (?population > 9999)
  # minus { ?place wdt:P31 wd:Q3024240 } .
  ?place rdfs:label ?placename .
  FILTER (lang(?placename)='en')
} 
"
