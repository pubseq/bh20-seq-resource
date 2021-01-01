import os

def is_integer(string_to_check):
    try:
        int(string_to_check)
        return True
    except ValueError:
        return False

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def check_and_get_ontology_dictionaries(dir_ontology_dictionaries):
    # Check duplicated entry looking at all dictionaries
    field_to_term_to_uri_dict = {}

    path_dict_xxx_csv_list = [os.path.join(dir_ontology_dictionaries, name_xxx_csv) for name_xxx_csv in
                              os.listdir(dir_ontology_dictionaries) if name_xxx_csv.endswith('.csv')]

    for path_dict_xxx_csv in path_dict_xxx_csv_list:
        print('Read {}'.format(path_dict_xxx_csv))

        with open(path_dict_xxx_csv) as f:
            for line in f:
                if len(line.split(',')) > 2:
                    term, uri = line.strip('\n').split('",')
                else:
                    term, uri = line.strip('\n').split(',')

                term = term.strip('"')

                if term in field_to_term_to_uri_dict:
                    print('Warning: in the dictionaries there are more entries for the same term ({}).'.format(term))
                    continue

                field_to_term_to_uri_dict[term] = uri

    # Prepare separated dictionaries (to avoid, for example, that a valid IRI for species is accepted as specimen)
    field_to_term_to_uri_dict = {}

    for path_dict_xxx_csv in path_dict_xxx_csv_list:
        field = os.path.basename(path_dict_xxx_csv).split('.')[0]

        field_to_term_to_uri_dict[field] = {}

        with open(path_dict_xxx_csv) as f:
            for line in f:
                if len(line.split(',')) > 2:
                    term, uri = line.strip('\n').split('",')
                else:
                    term, uri = line.strip('\n').split(',')

                term = term.strip('"')

                if term in field_to_term_to_uri_dict[field]:
                    print('Warning: in the {} dictionary there are more entries for the same term ({}).'.format(field, term))
                    continue

                field_to_term_to_uri_dict[field][term] = uri

    return field_to_term_to_uri_dict