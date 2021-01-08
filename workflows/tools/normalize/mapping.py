# Normalization steps
#
# This library contains generic logic to normalize (string) data and
# transforms strings to URIs.  It should be applicable to data from
# any source (GenBank, ENA etc).
#
#   Important: missing data should be missing or None! Do not fill
#   in data by 'guessing'.
#
#   When data is malformed a warning should be logged and added to the
#   warning list. Functions should be small enough to return only 1
#   warning!
#
#   Pjotr Prins (c) 2021

import re
import types

def host_species(host,mapping):
    Homo_sapiens = "http://purl.obolibrary.org/obo/NCBITaxon_9606"

    SPECIES_TERMS = { # since Python 3.7 dict is ordered! Note that re is allowed
        "human": Homo_sapiens,
        "sapiens": Homo_sapiens,
        "Mustela lutreola": "http://purl.obolibrary.org/obo/NCBITaxon_9666",
        "Manis javanica": "http://purl.obolibrary.org/obo/NCBITaxon_9974",
        "Felis catus": "http://purl.obolibrary.org/obo/NCBITaxon_9685",
        "Panthera tigris": "http://purl.obolibrary.org/obo/NCBITaxon_419130",
        "Canis lupus": "http://purl.obolibrary.org/obo/NCBITaxon_9615",
        # Mink:
        "vison": "http://purl.obolibrary.org/obo/NCBITaxon_452646"
        }

    warning = None
    host = types.SimpleNamespace(**host)
    if not 'obolibrary' in host.host_species:
        key = host.host_species
        host.host_species = None
        if key in mapping:
            host.host_species = mapping[key]
        else:
            for term in SPECIES_TERMS:
                p = re.compile(".*?"+term,re.IGNORECASE)
                m = p.match(key)
                if m: host.host_species = SPECIES_TERMS[term]
        if not host.host_species:
            warning = f"No URI mapping for host_species <{key}>"
        if host.host_species == Unknown or host.host_species == None:
            del(host.host_species)
    return host.__dict__,warning

Unknown = "Not found" # So as not to create a warning

def specimen_source(sample,mapping):
    Oronasopharynx = "http://purl.obolibrary.org/obo/NCIT_C155835"
    Oropharyngeal = "http://purl.obolibrary.org/obo/NCIT_C155835"
    Nasopharyngeal = "http://purl.obolibrary.org/obo/NCIT_C155831"
    Bronchoalveolar_Lavage_Fluid = "http://purl.obolibrary.org/obo/NCIT_C13195"
    Saliva = "http://purl.obolibrary.org/obo/NCIT_C13275"
    Nasal_Swab = Nasopharyngeal # "http://purl.obolibrary.org/obo/NCIT_C132119"
    Frozen_Food = "https://www.wikidata.org/wiki/Q751728"
    Bronchoalveolar_Lavage = "http://purl.obolibrary.org/obo/NCIT_C13195",
    Biospecimen = "http://purl.obolibrary.org/obo/NCIT_C70699"
    SPECIMEN_TERMS = { # since Python 3.7 dict is ordered! Note that re is allowed
        "Oronasopharynx": Oronasopharynx,
        "orophar": Oropharyngeal,
        "pharyngeal": Nasopharyngeal,
        "\snares": Nasal_Swab,
        "saliva": Saliva,
        "swab": Nasal_Swab,
        "broncho": Bronchoalveolar_Lavage,
        "seafood": Frozen_Food,
        "packaging": Frozen_Food,
        "specimen": Biospecimen,
        "patient": Biospecimen,
        "uknown": Unknown,
        "unknown": Unknown
        }
    warning = None
    sample = types.SimpleNamespace(**sample)
    try:
        if sample.specimen_source:
            keys = sample.specimen_source
            sample.specimen_source = []
            for key in keys:
                if 'obolibrary' in key:
                    sample.specimen_source.append(key)
                    continue
                if key in mapping:
                    sample.specimen_source.append(mapping[key])
                else:
                    for term in SPECIMEN_TERMS:
                        p = re.compile(".*?"+term,re.IGNORECASE)
                        m = p.match(key)
                        if m: sample.specimen_source = [SPECIMEN_TERMS[term]]
                if len(sample.specimen_source)==0:
                    warning = f"No URI mapping for specimen_source <{key}>"
        if sample.specimen_source == Unknown or sample.specimen_source == None:
            del(sample.specimen_source)
    except AttributeError:
        pass
    return sample.__dict__,warning
