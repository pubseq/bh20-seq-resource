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
    warning = None
    host = types.SimpleNamespace(**host)
    if not 'obolibrary' in host.host_species:
        key = host.host_species
        if key in mapping:
            host.host_species = mapping[key]
        else:
            warning = f"No URI mapping for host_species <{key}>"
    return host.__dict__,warning

Unknown = "Not found" # So as not to create a warning

def specimen_source(sample,mapping):
    Oronasopharynx = "http://purl.obolibrary.org/obo/NCIT_C155835"
    Oropharyngeal = "http://purl.obolibrary.org/obo/NCIT_C155835"
    Nasopharyngeal = "http://purl.obolibrary.org/obo/NCIT_C155831"
    Bronchoalveolar_Lavage_Fluid = "http://purl.obolibrary.org/obo/NCIT_C13195"
    Saliva = "http://purl.obolibrary.org/obo/NCIT_C13275"
    Nasal_Swab = "http://purl.obolibrary.org/obo/NCIT_C132119"
    Frozen_Food = "https://www.wikidata.org/wiki/Q751728"
    SPECIMEN_TERMS = { # since Python 3.7 dict is ordered! Note that re is allowed
        "Oronasopharynx": Oronasopharynx,
        "orophar": Oropharyngeal,
        "pharyngeal": Nasopharyngeal,
        "\snares": Nasal_Swab,
        "saliva": Saliva,
        "swab": Nasal_Swab,
        "seafood": Frozen_Food,
        "packaging": Frozen_Food,
        "uknown": Unknown,
        "unknown": Unknown
        }
    warning = None
    sample = types.SimpleNamespace(**sample)
    try:
        if sample.specimen_source and \
           not 'obolibrary' in sample.specimen_source and \
           not 'wikidata' in sample.specimen_source:
            key = sample.specimen_source
            sample.specimen_source = None
            if key in mapping:
                sample.specimen_source = mapping[key]
            else:
                for term in SPECIMEN_TERMS:
                    p = re.compile(".*?"+term,re.IGNORECASE)
                    m = p.match(key)
                    if m: sample.specimen_source = SPECIMEN_TERMS[term]
        if not sample.specimen_source:
            warning = f"No URI mapping for specimen_source <{key}>"
        if sample.specimen_source == Unknown or sample.specimen_source == None:
            del(sample.specimen_source)
    except AttributeError:
        pass
    return sample.__dict__,warning
