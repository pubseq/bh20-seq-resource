# ---- BELOW IS JUST FOR REFERENCE ----

if field_in_yaml == 'sequencing_coverage':
    # A regular expression would be better!
    try:
        technology[field_in_yaml] = [
            float(tech_info_to_parse.replace('(average)', '').replace("reads/nt", '').
                  replace('(average for 6 sequences)', '').replace(',', '.').strip(' xX>'))
        ]
    except ValueError:
        print(accession_version, "Couldn't make sense of Coverage '%s'" % tech_info_to_parse)
        pass


elif GBQualifier_name_text == 'collected_by':
    if any([x in GBQualifier_value_text.lower() for x in ['institute', 'hospital', 'city', 'center']]):
        sample['collecting_institution'] = GBQualifier_value_text
    else:
        sample['collector_name'] = GBQualifier_value_text

elif GBQualifier_name_text == 'isolation_source':
if GBQualifier_value_text.upper() in field_to_term_to_uri_dict['ncbi_speciesman_source']:
    GBQualifier_value_text = GBQualifier_value_text.upper()  # For example, in case of 'usa: wa'
