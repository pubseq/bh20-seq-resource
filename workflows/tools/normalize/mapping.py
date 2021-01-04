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

def specimen_source(sample,mapping):
    warning = None
    sample = types.SimpleNamespace(**sample)
    try:
        if sample.specimen_source and not 'obolibrary' in sample.specimen_source:
            key = sample.specimen_source
            if key in mapping:
                sample.specimen_source = mapping[key]
            else:
                sample.specimen_source = None
                warning = f"No URI mapping for specimen_source <{key}>"
    except AttributeError:
        pass
    if not sample.specimen_source: del(sample.specimen_source)
    return sample.__dict__,warning
