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

# Little cleaning
GBQualifier_value_text = GBQualifier_value_text.strip("/'")

if GBQualifier_value_text in field_to_term_to_uri_dict['ncbi_speciesman_source']:
    sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source'][GBQualifier_value_text]]
else:
    if GBQualifier_value_text.lower() in ['np/op', 'np-op', 'np/op swab', 'np/np swab', 'nasopharyngeal and oropharyngeal swab', 'nasopharyngeal/oropharyngeal swab', 'combined nasopharyngeal and oropharyngeal swab', 'naso and/or oropharyngeal swab']:
        sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasopharyngeal swab'], field_to_term_to_uri_dict['ncbi_speciesman_source']['oropharyngeal swab']]
    elif GBQualifier_value_text.lower() in ['nasopharyngeal swab/throat swab', 'nasopharyngeal/throat swab', 'nasopharyngeal swab and throat swab', 'nasal swab and throat swab', 'nasopharyngeal aspirate/throat swab', 'Nasopharyngeal/Throat']:
        sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasopharyngeal swab'], field_to_term_to_uri_dict['ncbi_speciesman_source']['throat swab']]
    elif GBQualifier_value_text.lower() in ['nasopharyngeal aspirate & throat swab', 'nasopharyngeal aspirate and throat swab']:
        sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasopharyngeal aspirate'], field_to_term_to_uri_dict['ncbi_speciesman_source']['throat swab']]
    elif GBQualifier_value_text.lower() in ['nasal swab and throat swab']:
        sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasal swab'], field_to_term_to_uri_dict['ncbi_speciesman_source']['throat swab']]
    elif GBQualifier_value_text.lower() in ['nasal-swab and oro-pharyngeal swab']:
        sample['specimen_source'] = [field_to_term_to_uri_dict['ncbi_speciesman_source']['nasal swab'], field_to_term_to_uri_dict['ncbi_speciesman_source']['oropharyngeal swab']]
    else:
        missing_value_list.append('\t'.join([accession_version, 'specimen_source', GBQualifier_value_text]))
