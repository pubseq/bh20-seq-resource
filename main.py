import tempfile
import shutil
import subprocess
import os
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

@app.route('/')
def send_form():
    """
    Send the file upload form/front page.
    """
    return send_from_directory('pages', 'index.html')
    
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
    
    
@app.route('/submit', methods=['POST'])
def recieve_files():
    """
    Recieve the uploaded files.
    """
    
    # We're going to work in one directory per request
    dest_dir = tempfile.mkdtemp()
    try:
    
        print(request)
        print(request.files)
    
        if 'fasta' not in request.files:
            return (render_template('error.html',
                error_message="You did not include a FASTA file."), 403)
        if 'metadata' not in request.files:
            return (render_template('error.html',
                error_message="You did not include a metadata file."), 403)
        
        fasta_dest = os.path.join(dest_dir, 'fasta.fa')
        metadata_dest = os.path.join(dest_dir, 'metadata.json')
                
        try:
            with open(fasta_dest, 'wb') as out_stream:
                copy_with_limit(request.files.get('fasta').stream, out_stream)
            with open(metadata_dest, 'wb') as out_stream:
                copy_with_limit(request.files.get('metadata').stream, out_stream)
        except FileTooBigError as e:
            # Delegate to the 413 error handler
            return handle_large_file(e)
            
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
        
    
    

