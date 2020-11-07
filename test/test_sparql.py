# Run with python3 test/test_sparql.py

import unittest
import requests
import logging



class TestSPARQL(unittest.TestCase):

    def test_sparql(self):
        sparqlURL='http://sparql.genenetwork.org/sparql/'
        id = "http://collections.lugli.arvadosapi.com/c=0002e93b86ad77824620bf938b97e134+126/sequence.fasta"
        id = "MT800005.1"
        query=f"""
PREFIX pubseq: <http://biohackathon.org/bh20-seq-schema#MainSchema/>
PREFIX sio: <http://semanticscience.org/resource/>
select distinct ?sample ?geoname ?date ?source ?geo ?sampletype ?institute ?sequenceuri
{{
   ?sample sio:SIO_000115 "{id}" .
   ?sequenceuri pubseq:sample ?sample .
   ?sample <http://purl.obolibrary.org/obo/GAZ_00000448> ?geo .
   ?geo rdfs:label ?geoname .
   ?sample <http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#C25164> ?date .
   OPTIONAL {{ ?sample <http://edamontology.org/data_2091> ?source }}
   OPTIONAL {{ ?sample <http://purl.obolibrary.org/obo/OBI_0001479> ?sampletype }}
   OPTIONAL {{ ?sample <http://purl.obolibrary.org/obo/NCIT_C41206> ?institute }}
}}
        """
        print(query)
        payload = {'query': query, 'format': 'json'}
        r = requests.get(sparqlURL, params=payload)
        result = r.json()['results']['bindings']
        # for now we just take the first one
        print(result)
        self.assertEqual(result[0]['geoname']['value'],'Mahuva')

if __name__ == '__main__':
    unittest.main()
