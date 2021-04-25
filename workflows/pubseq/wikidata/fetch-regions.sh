# curl -G https://query.wikidata.org/sparql -H "Accept: text/tab-separated-values; charset=utf-8" --data-urlencode query="
curl -G https://query.wikidata.org/sparql -H "Accept: text/csv; charset=utf-8" --data-urlencode query="
select distinct ?placename ?place ?country ?population ?coor where {
  VALUES ?v { wd:Q82794 wd:Q107390 wd:Q34876 wd:Q9316670 wd:Q515 }
  ?statetype wdt:P279+ ?v .
  ?place wdt:P31 ?statetype ;
         wdt:P17 ?country ;
         wdt:P625 ?coor;
         wdt:P1082 ?population .
  FILTER (?population > 99999)
  ?place rdfs:label ?placename .
  FILTER (lang(?placename)='en')
}

"
