# This script by @LLTommy queries the main SPARQL end point to find what
# collections are missing country information for GPS coordinates, such
#
# <http://www.wikidata.org/entity/Q657004> rdfs:label "Canterbury Region" ;
#    ns1:P17 <http://www.wikidata.org/entity/Q664> ;
#    ns1:P625 "Point(172.0 -43.6)" .
#
# See also the ./readme.md

import requests
import csv
from rdflib import Graph, Literal, RDF, URIRef
import time


sparqlHeaders={'User-Agent': 'genenetworkCrawler/1.0 (covid-19.genenetwork.org; pjotr.public821@thebird.nl) genenetworkCrawler/1.0'}


def callSPARQL(query):
    payload = {'query': query, 'format': 'json'}
    try:
        r = requests.get('https://query.wikidata.org/sparql', params=payload, headers=sparqlHeaders)
        time.sleep(1)
        # Slow process down, in case we did send too many processes. Sleep, then try again
        if 'Retry_After' in r.headers:
            print(r.headers)
            time.sleep(45)
            r = requests.get('https://query.wikidata.org/sparql', params=payload, headers=sparqlHeaders)

        result = r.json()['results']['bindings']
    except Exception as e:
        print("Error during SPARQL call. We abort the process and have to investigate")
        print(r)
        print(r.content)
        print(r.url)
        raise Exception(e)
    return result


g = Graph()

test_query="""
# Use with https://query.wikidata.org/
SELECT DISTINCT ?a ?label ?country ?continent ?coor WHERE {
    BIND (XXX as ?a) .
    OPTIONAL {
        ?a wdt:P625 ?coor.
    }
    ?a rdfs:label ?label .
    ?a wdt:P17 ?country.
    ?country rdfs:label ?country_label .
    OPTIONAL {
        ?country wdt:P30 ?continent.
        ?continent rdfs:label ?continent_label
        FILTER (lang(?continent_label)='en')
    }
    FILTER (lang(?country_label)='en')
    FILTER (lang(?label)='en')
}
"""

# wdt:P625 are GEO coordinates

query = """
construct {
    ?a wdt:P625 ?c.
    ?a rdfs:label ?label .
    ?a wdt:P17 ?country.
    ?country rdfs:label ?country_label .
    ?country wdt:P30 ?continent .
    ?continent rdfs:label ?continent_label .
} WHERE
{
    BIND (XXX as ?a) .
    ?a rdfs:label ?label .
    ?a wdt:P17 ?country.
    ?country rdfs:label ?country_label .
    ?country wdt:P30 ?continent.
    ?continent rdfs:label ?continent_label
    FILTER (lang(?continent_label)='en')
    FILTER (lang(?country_label)='en')
    FILTER (lang(?label)='en')
    OPTIONAL { ?a wdt:P625 ?c }
}
"""""

outputFile = 'input_location.csv'

with open(outputFile, 'r') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
    counter=0
    for row in spamreader:
        counter=counter+1

        try:
            testq = test_query.replace("XXX", "<"+row[0]+">")
            print(testq)
            tmpquery=query.replace("XXX", "<"+row[0]+">")
            print(tmpquery)

            returnReply=callSPARQL(tmpquery)
            print(returnReply)

            for row in returnReply:
                print(row)
                sub=URIRef(row['subject']['value'])
                pred=URIRef(row['predicate']['value'])

                if row['object']['type'] == 'uri':
                    obj =  URIRef(row['object']['value'])

                elif row['object']['type'] == 'literal':
                    obj= Literal(row['object']['value'])
                g.add(( sub, pred, obj ))

        except Exception as e:
            print(e)
            raise

print(g.serialize(format='n3').decode("utf-8"))
g.serialize(destination='enriched_output.txt', format='turtle')
