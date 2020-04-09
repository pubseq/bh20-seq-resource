import collections
import tempfile
import shutil
import subprocess
import os
import re
import string
import yaml
import urllib.request
from flask import Flask, request, redirect, send_file, send_from_directory, render_template

app = Flask(__name__, static_url_path='/static', static_folder='static')

# Limit file upload size. We shouldn't be working with anything over 1 MB; these are small genomes.
# We will enforce the limit ourselves and set a higher safety limit here.
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024

# When a file is too big we get a 413.
@app.errorhandler(413)
def handle_large_file(e):
    return (render_template('error.html',
        error_message="One of your files is too large. The maximum file size is 1 megabyte."), 413)


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

    return string.capwords(field_name.replace('_', ' '))

def generate_form(schema):
    """
    Linearize the schema and send a bunch of dicts.
    Each dict either has a 'heading' (in which case we put a heading for a
    form section in the template) or an 'id', 'label', 'type', and 'required'
    (in which case we make a form field in the template).
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

            ref_url = None
            if not isinstance(field_type, str):
                # If the type isn't a string
                # See if it has a more info/what goes here URL
                ref_url = field_type.get('jsonldPredicate', {}).get('_id', None)
                # Grab out its type field
                field_type = field_type.get('type', '')

            # Decide if the field is optional (type ends in ?)
            optional = False
            if len(field_type) > 0 and field_type[-1] == '?':
                # It's optional
                optional = True
                # Drop the ?
                field_type = field_type[:-1]

            if field_type in by_name:
                # This is a subrecord. We need to recurse
                for item in walk_fields(field_type, parent_keys + [field_name], subtree_optional or optional):
                    yield item
            else:
                # We know how to make a string input
                record = {}
                record['id'] = '.'.join(parent_keys + [field_name])
                record['label'] = name_to_label(field_name)
                record['required'] = not optional and not subtree_optional
                if ref_url:
                    record['ref_url'] = ref_url
                if field_type == 'string':
                    record['type'] = 'text' # HTML input type
                elif field_type == 'int':
                    record['type'] = 'number'
                else:
                    raise NotImplementedError('Unimplemented field type {} in {} in metadata schema'.format(field_type, type_name))
                yield record

    return list(walk_fields(root_name))

# At startup, we need to load the current metadata schema so we can make a form for it
METADATA_SCHEMA = yaml.safe_load(urllib.request.urlopen('https://raw.githubusercontent.com/arvados/bh20-seq-resource/master/bh20sequploader/bh20seq-schema.yml'))
FORM_ITEMS = generate_form(METADATA_SCHEMA)

@app.route('/')
def send_form():
    """
    Send the file upload form/front page.
    """

    return render_template('form.html', fields=FORM_ITEMS)

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

def parse_input(input_string, html_type):
    """
    Parse an input from the given HTML input type into a useful Python type.

    Raise ValueError if something does not parse.
    Raise NotImplementedError if we forgot to implement a type.
    """

    if html_type == 'text':
        return input_string
    elif html_type == 'number':
        return int(input_string)
    else:
        raise NotImplementedError('Unimplemented input type: {}'.format(html_type))

@app.route('/submit', methods=['POST'])
def receive_files():
    """
    Receive the uploaded files.
    """

    # We're going to work in one directory per request
    dest_dir = tempfile.mkdtemp()
    fasta_dest = os.path.join(dest_dir, 'fasta.fa')
    metadata_dest = os.path.join(dest_dir, 'metadata.json')
    try:
        if 'fasta' not in request.files:
            return (render_template('error.html',
                error_message="You did not include a FASTA file."), 403)
        try:
            with open(fasta_dest, 'wb') as out_stream:
                copy_with_limit(request.files.get('fasta').stream, out_stream)
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

            for item in FORM_ITEMS:
                # Pull all the field values we wanted from the form
                if 'heading' in item:
                    continue

                if item['id'] in request.form and len(request.form[item['id']]) > 0:
                    # We have this thing. Make a place in the dict tree for it.
                    parts = item['id'].split('.')
                    key = parts[-1]
                    # Remove leading 'metadata'
                    path = parts[1:-1]
                    dest_dict = metadata
                    for parent in path:
                        if parent not in dest_dict:
                            dest_dict[parent] = {}
                        dest_dict = dest_dict[parent]

                    try:
                        # Now finally add the item
                        dest_dict[key] = parse_input(request.form[item['id']], item['type'])
                    except ValueError:
                        # We don't like that input
                        return (render_template('error.html',
                            error_message="You provided an unacceptable value for the metadata item {}".format(item['id'])), 403)
                elif item['required']:
                    return (render_template('error.html',
                            error_message="You omitted the required metadata item {}".format(item['id'])), 403)

            # Now serialize the file with all the items
            with open(metadata_dest, 'w') as out_stream:
                yaml.dump(metadata, out_stream)
        else:
            return (render_template('error.html',
                    error_message="You did not include metadata."), 403)

        # Try and upload files to Arvados
        result = subprocess.run(['bh20-seq-uploader', fasta_dest, metadata_dest],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.returncode != 0:
            # It didn't work. Complain.
            error_message="Upload failed. Uploader returned {} and said:\n{}".format(result.returncode, result.stderr)
            return (render_template('error.html', error_message=error_message), 403)
        else:
            # It worked. Say so.
            return render_template('success.html', log=result.stdout.decode('utf-8', errors='replace'))
    finally:
        shutil.rmtree(dest_dir)
