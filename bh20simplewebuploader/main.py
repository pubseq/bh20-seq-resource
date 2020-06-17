import collections
import itertools
import tempfile
import shutil
import subprocess
import logging
import os
import sys
import re
import string
import yaml
import pkg_resources
from flask import Flask, request, redirect, send_file, send_from_directory, render_template, jsonify
import os.path
import requests

logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__ )
log.debug("Entering web uploader")

if not os.path.isfile('bh20sequploader/mainx.py'):
    print("WARNING: run FLASK from the root of the source repository!", file=sys.stderr)

app = Flask(__name__, static_url_path='/static', static_folder='static')

# Limit file upload size. We shouldn't be working with anything over 1 MB; these are small genomes.
# We will enforce the limit ourselves and set a higher safety limit here.
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024

# When a file is too big we get a 413.
@app.errorhandler(413)
def handle_large_file(e):
    return (render_template('error.html',
        error_message="One of your files is too large. The maximum file size is 50 megabytes."), 413)


def type_to_heading(type_name):
    """
    Turn a type name like "sampleSchema" from the metadata schema into a human-readable heading.
    """

    # Remove camel case
    decamel = re.sub('([A-Z])', r' \1', type_name)
    # Split
    parts = decamel.split()
    # Capitalize words and remove unwanted components
    filtered = [part.capitalize() for part in parts if (part.lower() != 'schema' and part != '')]
    # Reassemble
    return ' '.join(filtered)

def name_to_label(field_name):
    """
    Turn a filed name like "host_health_status" from the metadata schema into a human-readable label.
    """

    # May end in a number, which should be set off by a space
    set_off_number = re.sub('([0-9]+)$', r' \1', field_name)

    return string.capwords(set_off_number.replace('_', ' '))

def is_iri(string):
    """
    Return True if the given string looks like an IRI, and False otherwise.

    Used for finding type IRIs in the schema.

    Right now only supports http(s) URLs because that's all we have in our schema.
    """

    return string.startswith('http')

def generate_form(schema, options):
    """
    Linearize the schema into a list of dicts.

    Each dict either has a 'heading' (in which case we put a heading for a
    form section in the template) or an 'id', 'label', 'type', and 'required'
    (in which case we make a form field in the template).

    Non-heading dicts with type 'select' will have an 'options' field, with a
    list of (name, value) tuples, and represent a form dropdown element.

    Non-heading dicts with type 'number' may have a 'step', which, if <1 or
    'any', allows the number to be a float.

    Non-heading dicts may have a human-readable 'docstring' field describing
    them.

    Takes the deserialized metadata schema YAML, and also a deserialized YAML
    of option values. The option values are keyed on (unscoped) field name in
    the schema, and each is a dict of human readable option -> corresponding
    IRI.
    """

    # Get the list of form components, one of which is the root
    components = schema.get('$graph', [])

    # Find the root
    root_name = None
    # And also index components by type name
    by_name = {}
    for component in components:
        # Get the name of each
        component_name = component.get('name', None)
        if isinstance(component_name, str):
            # And remember how to map back form it
            by_name[component_name] = component
        if component.get('documentRoot', False):
            # Find whichever one is the root
            root_name = component_name


    def walk_fields(type_name, parent_keys=['metadata'], subtree_optional=False):
        """
        Do a traversal of the component tree.
        Yield a bunch of form item dicts, in order.
        Form IDs are .-separated keypaths for where they are in the structure.
        parent_keys is the path of field names to where we are in the root record's document tree.
        """

        if len(parent_keys) > 1:
            # First make a heading, if we aren't the very root of the form
            yield {'heading': type_to_heading(type_name)}

        for field_name, field_type in by_name.get(type_name, {}).get('fields', {}).items():
            # For each field

            ref_iri = None
            docstring = None
            if not isinstance(field_type, str):
                # If the type isn't a string

                # It may have documentation
                docstring = field_type.get('doc', None)

                # See if it has a more info/what goes here URL
                predicate = field_type.get('jsonldPredicate', {})
                # Predicate may be a URL, a dict with a URL in _id, maybe a
                # dict with a URL in _type, or a dict with _id and _type but no
                # URLs anywhere. Some of these may not technically be allowed
                # by the format, but if they occur, we might as well try to
                # handle them.
                if isinstance(predicate, str):
                    if is_iri(predicate):
                        ref_iri = predicate
                else:
                    # Assume it's a dict. Look at the fields we know about.
                    for field in ['_id', 'type']:
                        field_value = predicate.get(field, None)
                        if isinstance(field_value, str) and is_iri(field_value) and ref_iri is None:
                            # Take the first URL-looking thing we find
                            ref_iri = field_value
                            break


                # Now overwrite the field type with the actual type string
                field_type = field_type.get('type', '')

            # Decide if the field is optional (type ends in ?)
            optional = False
            if field_type.endswith('?'):
                # It's optional
                optional = True
                # Drop the ?
                field_type = field_type[:-1]

            # Decide if the field is a list (type ends in [])
            is_list = False
            if field_type.endswith('[]'):
                # It's a list
                is_list = True
                # Reduce to the normal type
                field_type = field_type[:-2]

            if field_type in by_name:
                # This is a subrecord. We need to recurse
                for item in walk_fields(field_type, parent_keys + [field_name], subtree_optional or optional):
                    yield item
            else:
                # This is a leaf field. We need an input for it.
                record = {}
                record['id'] = '.'.join(parent_keys + [field_name])
                record['label'] = name_to_label(field_name)
                record['required'] = not optional and not subtree_optional
                record['list'] = is_list
                if ref_iri:
                    record['ref_iri'] = ref_iri
                if docstring:
                    record['docstring'] = docstring

                if field_name in options:
                    # The field will be a 'select' type no matter what its real
                    # data type is.
                    record['type'] = 'select' # Not a real HTML input type. It's its own tag.
                    # We have a set of values to present
                    record['options'] = []
                    for name, value in options[field_name].items():
                        # Make a tuple for each one
                        record['options'].append((name, value))
                elif field_type == 'string':
                    if field_name.endswith('date'):
                        # Use a date picker to generate a good string.
                        # Comes back YYYY-MM-DD.
                        record['type'] = 'date'
                    else:
                        # Normal text string
                        record['type'] = 'text'
                elif field_type == 'int':
                    record['type'] = 'number'
                elif field_type == 'float' or field_type == 'double':
                    record['type'] = 'number'
                    # Choose a reasonable precision for the control
                    record['step'] = '0.0001'
                else:
                    raise NotImplementedError('Unimplemented field type {} in {} in metadata schema'.format(field_type, type_name))
                yield record

    return list(walk_fields(root_name))


# At startup, we need to load the metadata schema from the uploader module, so we can make a form for it
METADATA_SCHEMA = yaml.safe_load(pkg_resources.resource_stream("bh20sequploader", "bh20seq-schema.yml"))
METADATA_OPTION_DEFINITIONS = yaml.safe_load(pkg_resources.resource_stream("bh20sequploader", "bh20seq-options.yml"))
FORM_ITEMS = generate_form(METADATA_SCHEMA, METADATA_OPTION_DEFINITIONS)

@app.route('/')
def send_form():
    """
    Send the file upload form/front page.
    """

    return render_template('form.html', fields=FORM_ITEMS, menu='HOME')

class FileTooBigError(RuntimeError):
    """
    Raised when the user gives a file that is too large.
    """
    pass

def copy_with_limit(in_file, out_file, limit=1024*1024):
    """
    Copy a file stream, and raise FileTooBigError if the file is too big.
    """

    bytes_used = 0
    buf_size = 65536

    buf = in_file.read(buf_size)
    bytes_used += len(buf)
    while buf:
        if bytes_used > limit:
            raise FileTooBigError('Hit file length limit')
        out_file.write(buf)
        buf = in_file.read(buf_size)
        bytes_used += len(buf)

def parse_input(input_string, html_type, number_step=None):
    """
    Parse an input from the given HTML input type into a useful Python type.
    Also needs the step we sent to distinguish int fields and float/double fields.

    Raise ValueError if something does not parse.
    Raise NotImplementedError if we forgot to implement a type.
    """

    if html_type == 'text' or html_type == 'select':
        return input_string
    elif html_type == 'number':
        # May be an int or a float.
        if number_step is None:
            # TODO: Assumes we only use the step for floats
            return int(input_string)
        else:
            return float(input_string)
    elif html_type == 'date':
        # Don't do our own date validation; pass it on as a string
        return input_string
    else:
        raise NotImplementedError('Unimplemented input type: {}'.format(html_type))

@app.route('/submit', methods=['POST'])
def receive_files():
    """
    Receive the uploaded files.
    """

    # We're going to work in one directory per request
    dest_dir = tempfile.mkdtemp()
    # The uploader will happily accept a FASTQ with this name
    fasta_dest = os.path.join(dest_dir, 'fasta.fa')
    metadata_dest = os.path.join(dest_dir, 'metadata.json')
    try:
        if 'fasta' not in request.files:
            return (render_template('error.html',
                error_message="You did not include a FASTA or FASTQ file."), 403)
        try:
            with open(fasta_dest, 'wb') as out_stream:
                # Use a plausible file size limit for a little FASTQ
                copy_with_limit(request.files.get('fasta').stream, out_stream, limit=50*1024*1024)
        except FileTooBigError as e:
            # Delegate to the 413 error handler
            return handle_large_file(e)

        if request.form.get('metadata_type', None) == 'upload':
            if 'metadata' not in request.files:
                return (render_template('error.html',
                    error_message="You did not include a metadata file."), 403)
            try:
                with open(metadata_dest, 'wb') as out_stream:
                    copy_with_limit(request.files.get('metadata').stream, out_stream)
            except FileTooBigError as e:
                # Delegate to the 413 error handler
                return handle_large_file(e)
        elif request.form.get('metadata_type', None) == 'fill':
            # Build a metadata dict
            metadata = {}

            # When we have metadata for an item, use this to set it.
            # If it is an item in a list, set is_list=True
            def set_metadata(item_id, value, is_list=False):
                # We have this thing. Make a place in the dict tree for it.
                parts = item_id.split('.')
                key = parts[-1]
                # Remove leading 'metadata'
                path = parts[1:-1]
                dest_dict = metadata
                for parent in path:
                    if parent not in dest_dict:
                        dest_dict[parent] = {}
                    dest_dict = dest_dict[parent]

                if not is_list:
                    dest_dict[key] = value
                else:
                    if key not in dest_dict:
                        dest_dict[key] = []
                    dest_dict[key].append(value)

            for item in FORM_ITEMS:
                # Pull all the field values we wanted from the form
                if 'heading' in item:
                    continue

                if item['list']:
                    # This is a list, serialized into form fields

                    # We count how many values we got
                    value_count = 0

                    for index in itertools.count():
                        # Get [0] through [n], until something isn't there.
                        entry_id = '{}[{}]'.format(item['id'], index)

                        if index == 1000:
                            # Don't let them provide too much stuff.
                            return (render_template('error.html',
                                    error_message="You provided an extremely large number of values for the metadata item {}".format(item['id'])), 403)

                        if entry_id in request.form:
                            if len(request.form[entry_id]) > 0:
                                # Put an entry in the list
                                try:
                                    # Parse the item
                                    parsed = parse_input(request.form[entry_id], item['type'], item.get('step', None))
                                except ValueError:
                                    # We don't like that input
                                    return (render_template('error.html',
                                            error_message="You provided an unacceptable value for the metadata item {}".format(entry_id)), 403)
                                # Save it
                                set_metadata(item['id'], parsed, is_list=True)
                                value_count += 1
                            else:
                                # Empty items are silently skipped.
                                pass
                        else:
                            # We have run out of form fields for this list.
                            break

                        if item['required'] and value_count == 0:
                            # They forgot a required item. Maybe all entries were empty.
                            return (render_template('error.html',
                                    error_message="You omitted any values for the required metadata item {}".format(item['id'])), 403)

                elif item['id'] in request.form and len(request.form[item['id']]) > 0:
                    # Not a list, but a single item which is present.
                    try:
                        # Parse the item
                        parsed = parse_input(request.form[item['id']], item['type'], item.get('step', None))
                    except ValueError:
                        # We don't like that input
                        return (render_template('error.html',
                            error_message="You provided an unacceptable value for the metadata item {}".format(item['id'])), 403)
                    # Save it
                    set_metadata(item['id'], parsed)
                elif item['required']:
                    return (render_template('error.html',
                            error_message="You omitted the required metadata item {}".format(item['id'])), 403)

            # Now serialize the file with all the items
            with open(metadata_dest, 'w') as out_stream:
                yaml.dump(metadata, out_stream)
        else:
            return (render_template('error.html',
                    error_message="You did not include metadata."), 403)

        # Try and upload files to Arvados using the sequence uploader CLI

        cmd = ['python3','bh20sequploader/main.py', fasta_dest, metadata_dest]
        print(" ".join(cmd),file=sys.stderr)
        result = subprocess.run(cmd,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.returncode != 0:
            # It didn't work. Complain.
            error_message="Uploader returned value {} and said:".format(result.returncode) + str(result.stderr.decode('utf-8'))
            print(error_message, file=sys.stderr)
            return (render_template('error.html', error_message=error_message), 403)
        else:
            # It worked. Say so.
            return render_template('success.html', log=result.stdout.decode('utf-8', errors='replace'))
    finally:
        shutil.rmtree(dest_dir)

def get_html_body(fn):
    buf = ""
    in_body = False
    begin_body = re.compile(r"<body>",re.IGNORECASE)
    end_body = re.compile(r"(</body>|.*=\"postamble\")",re.IGNORECASE)
    with open(fn) as f:
        for line in f:
            if end_body.match(line):
                break
            if in_body:
                buf += line
            elif begin_body.match(line):
                in_body = True
    return buf

@app.route('/download')
def download_page():
    buf = get_html_body('doc/web/download.html')
    return render_template('about.html',menu='DOWNLOAD',embed=buf)

@app.route('/demo')
def demo_page():
    return render_template('demo.html',menu='DEMO')

@app.route('/blog',methods=['GET'])
def blog_page():
    blog_content = request.args.get('id') # e.g. using-covid-19-pubseq-part3
    buf = None;
    if blog_content:
        buf = get_html_body('doc/blog/'+blog_content+'.html')
    return render_template('blog.html',menu='BLOG',embed=buf,blog=blog_content)


@app.route('/about')
def about_page():
    buf = get_html_body('doc/web/about.html')
    return render_template('about.html',menu='ABOUT',embed=buf)

##
@app.route('/map')
def map_page():
    return render_template('map.html',menu='DEMO')



## Dynamic API functions starting here
## This is quick and dirty for now, just to get something out and demonstrate the queries
## Feel free to rename the functions/endpoints, feel free to process result so we get nicer JSON
## but most likley you don't want to touch the queries, Cheers.
baseURL='http://sparql.genenetwork.org/sparql/'

@app.route('/api/getCount', methods=['GET'])
def getCount():
    query="""
PREFIX pubseq: <http://biohackathon.org/bh20-seq-schema#MainSchema/>
select (COUNT(distinct ?dataset) as ?num)
{
   ?dataset pubseq:submitter ?id .
   ?id ?p ?submitter
}
"""
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    # [{'num': {'type': 'typed-literal', 'datatype': 'http://www.w3.org/2001/XMLSchema#integer', 'value': '1352'}}]
    # print(result, file=sys.stderr)
    return jsonify({'sequences': int(result[0]["num"]["value"])})

@app.route('/api/getAllaccessions', methods=['GET'])
def getAllaccessions():
    query="""SELECT DISTINCT ?fasta ?value WHERE {?fasta ?x[ <http://edamontology.org/data_2091> ?value ]}"""
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return jsonify([{'uri': x['fasta']['value'],
                     'value': x['value']['value']} for x in result])


# parameter must be encoded e.g. http://arvados.org/keep:6e6276698ed8b0e6cd21f523e4f91179+123/sequence.fasta must become
# http%3A%2F%2Fcollections.lugli.arvadosapi.com%2Fc%3D00a6af865453564f6a59b3d2c81cc7c1%2B123%2Fsequence.fasta
@app.route('/api/getDetailsForSeq', methods=['GET'])
def getDetailsForSeq():
    seq_id = request.args.get('seq')
    query="""SELECT DISTINCT ?key ?key_label ?value WHERE {
    <placeholder> ?x [?key ?value] .
    OPTIONAL {?key <http://www.w3.org/2000/01/rdf-schema#label> ?key_tmp_label } .
    BIND(IF(BOUND(?key_tmp_label),?key_tmp_label, ?key) as ?key_label)}"""
    query=query.replace("placeholder", seq_id)
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return jsonify([{'uri': x['key']['value'], 'key_label': x['key_label']['value'],
                     'value': x['value']['value']} for x in result])

# Endpoint should provide all necessary information to draw a map (!)
@app.route('/api/getCountByGPS', methods=['GET'])
def getCountByGPS():
    query="""SELECT DISTINCT ?location ?location_label ?GPS (count(?fasta) as ?fastaCount) WHERE {
    ?fasta ?x[ <http://purl.obolibrary.org/obo/GAZ_00000448> ?location] .
    ?location <http://www.wikidata.org/prop/direct/P625> ?GPS .
    OPTIONAL { ?location rdfs:label ?key_tmp_label }
    BIND(IF(BOUND(?key_tmp_label),?key_tmp_label, ?location) as ?location_label)
    }
    GROUP BY ?location ?location_label ?GPS
    """
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return jsonify([{'count': x['fastaCount']['value'],
                     'Location': x['location']['value'],
                     'LocationLabel': x['location_label']['value'],
                     'GPS' :x['GPS']['value'][6:-1]} for x in result])


@app.route('/api/getSEQCountbytech', methods=['GET'])
def getSEQCountbytech():
    query="""SELECT ?tech ?tech_label (count(?fasta) as ?fastaCount) WHERE
    {?fasta ?x [<http://purl.obolibrary.org/obo/OBI_0600047>  ?tech] .
     OPTIONAL {?tech <http://www.w3.org/2000/01/rdf-schema#label> ?tech_tmp_label } .
     BIND(IF(BOUND(?tech_tmp_label), ?tech_tmp_label,?tech) as ?tech_label)}
    GROUP BY ?tech ?tech_label ORDER BY DESC (?fastaCount)
    """
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return jsonify([{'count': x['fastaCount']['value'],
                     'key': x['tech']['value'],
                     'label': x['tech_label']['value']} for x in result])

## List all Sequences/submissions by a given tech, as example e.g. http://purl.obolibrary.org/obo/OBI_0000759
## Has to be encoded again so should be --> http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FOBI_0000759
@app.route('/api/getSEQbytech', methods=['GET'])
def getSEQbytech():
    query="""SELECT ?fasta WHERE
    {?fasta ?x [<http://purl.obolibrary.org/obo/OBI_0600047>  <placeholder>] }
    """
    tech = request.args.get('tech')
    query=query.replace("placeholder", tech)
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return str(result)


## Example location, encoded http%3A%2F%2Fwww.wikidata.org%2Fentity%2FQ1223
@app.route('/api/getSEQbyLocation', methods=['GET'])
def getSEQbyLocation():
    query="""SELECT ?fasta WHERE {?fasta ?x[ <http://purl.obolibrary.org/obo/GAZ_00000448> <placeholder>]}"""
    location=request.args.get('location')
    query=query.replace("placeholder", location)
    print(query)
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return str(result)


@app.route('/api/getSEQCountbyLocation', methods=['GET'])
def getSEQCountbyLocation():
    query="""SELECT ?geoLocation ?geoLocation_label (count(?fasta) as ?fastaCount)  WHERE
    {?fasta ?x [<http://purl.obolibrary.org/obo/GAZ_00000448> ?geoLocation] .
    Optional {?geoLocation <http://www.w3.org/2000/01/rdf-schema#label> ?geoLocation_tmp_label}
    BIND(IF(BOUND(?geoLocation_tmp_label), ?geoLocation_tmp_label, ?geoLocation) as ?geoLocation_label)}
    GROUP BY ?geoLocation ?geoLocation_label ORDER BY DESC (?fastaCount)
    """
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return jsonify([{'count': x['fastaCount']['value'],
                     'key': x['geoLocation']['value'],
                     'label': x['geoLocation_label']['value']} for x in result])


@app.route('/api/getSEQCountbyContinent', methods=['GET'])
def getSEQCountbyContinent():
    query="""SELECT DISTINCT ?continent ?continent_label (count(?fasta) as ?fastaCount) WHERE {
    ?fasta ?x[ <http://purl.obolibrary.org/obo/GAZ_00000448> ?location] .
   ?location <http://www.wikidata.org/prop/direct/P17> ?country .
    ?country <http://www.wikidata.org/prop/direct/P30> ?continent .
    OPTIONAL { ?continent rdfs:label ?key_tmp_label }
    BIND(IF(BOUND(?key_tmp_label),?key_tmp_label, ?location) as ?continent_label)
    }
    GROUP BY ?continent ?continent_label
    """
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return jsonify([{'count': x['fastaCount']['value'],
                     'key': x['continent']['value'],
                     'label': x['continent_label']['value']} for x in result])



@app.route('/api/getSEQCountbyCountryContinent', methods=['GET'])
def getSEQCountbyCountryContinent():
    query="""SELECT DISTINCT ?location ?location_label (count(?fasta) as ?fastaCount) WHERE {
    ?fasta ?x[ <http://purl.obolibrary.org/obo/GAZ_00000448> ?location] .
    ?location <http://www.wikidata.org/prop/direct/P30> <placeholder> .
    OPTIONAL { ?location rdfs:label ?key_tmp_label }
    BIND(IF(BOUND(?key_tmp_label),?key_tmp_label, ?location) as ?location_label)
    }
    GROUP BY ?location ?location_label
    """
    continent = request.args.get('continent')
    query = query.replace("placeholder", continent)
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return jsonify([{'count': x['fastaCount']['value'],
                     'key': x['location']['value'],
                     'label': x['location_label']['value']} for x in result])



@app.route('/api/getSEQCountbySpecimenSource', methods=['GET'])
def getSEQCountbySpecimenSource():
    query="""SELECT ?specimen_source ?specimen_source_label (count(?fasta) as ?fastaCount)  WHERE
    {?fasta ?x [<http://purl.obolibrary.org/obo/OBI_0001479>  ?specimen_source]
     Optional { ?specimen_source <http://www.w3.org/2000/01/rdf-schema#label> ?specimen_source_tmp_label}
     BIND(IF(BOUND(?specimen_source_tmp_label), ?specimen_source_tmp_label ,?specimen_source) as ?specimen_source_label)}
     GROUP BY ?specimen_source ?specimen_source_label
     ORDER BY DESC (?fastaCount)
    """
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return jsonify([{'count': x['fastaCount']['value'],
                     'key': x['specimen_source']['value'],
                     'label': x['specimen_source_label']['value']} for x in result])

# Example specimen http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCIT_C155831
@app.route('/api/getSEQbySpecimenSource', methods=['GET'])
def getSEQBySpecimenSource():
    query="""SELECT ?fasta ?specimen_source ?specimen_source_label  WHERE
    {?fasta ?x [<http://purl.obolibrary.org/obo/OBI_0001479> <placeholder>]
    BIND (concat(?specimen_source,"_label") as ?specimen_source_label)}
    """
    specimen=request.args.get('specimen')
    query = query.replace("placeholder", specimen)
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return str(result)

#No data for this atm
@app.route('/api/getSEQCountbyHostHealthStatus', methods=['GET'])
def getSEQCountbyHostHealthStatus():
    query="""SELECT ?health_status ?health_status_label (count(?fasta) as ?fastaCount)  WHERE
    {?fasta ?x [<http://purl.obolibrary.org/obo/NCIT_C25688> ?health_status]
    BIND (concat(?health_status,"_label") as ?health_status_label)}
    GROUP BY ?health_status ?health_status_label
    ORDER BY DESC (?fastaCount)
    """
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return str(result)

@app.route('/api/getSEQbyLocationAndTech', methods=['GET'])
def getSEQbyLocationAndTech():
    query="""SELECT ?fasta WHERE { ?fasta ?x [
        <http://purl.obolibrary.org/obo/GAZ_00000448> <placeholderLoc>; <http://purl.obolibrary.org/obo/OBI_0600047>  <placeholderTech> ]}"""
    location=request.args.get('location')
    tech=request.args.get('tech')
    query=query.replace("placeholderLoc", location)
    query = query.replace("placeholderTech", tech)
    print(query)
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return str(result)


# Example Location http%3A%2F%2Fwww.wikidata.org%2Fentity%2FQ1223
# Example specimen http%3A%2F%2Fpurl.obolibrary.org%2Fobo%2FNCIT_C155831
@app.route('/api/getSEQbyLocationAndSpecimenSource', methods=['GET'])
def getSEQbyLocationAndSpecimenSource():
    query="""SELECT ?fasta WHERE { ?fasta ?x [
        <http://purl.obolibrary.org/obo/GAZ_00000448> <placeholderLoc>; <http://purl.obolibrary.org/obo/OBI_0001479>  <placeholderSpecimen> ]}
    """
    location = request.args.get('location')
    specimen = request.args.get('specimen')
    query = query.replace("placeholderLoc", location)
    query = query.replace("placeholderSpecimen", specimen)
    print(query)
    payload = {'query': query, 'format': 'json'}
    r = requests.get(baseURL, params=payload)
    result = r.json()['results']['bindings']
    return str(result)
