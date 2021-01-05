import gzip
import tempfile
import pkg_resources
import magic
import io

import sys

path_fasta = sys.argv[1]
format_to_check = sys.argv[2]
path_valid_formats = sys.argv[3]


# ../../bh20sequploader/validation/formats

schema_resource = pkg_resources.resource_stream(__name__, path_valid_formats)
with tempfile.NamedTemporaryFile() as tmp:
    tmp.write(schema_resource.read())
    tmp.flush()
    check_format = magic.Magic(magic_file=tmp.name, uncompress=False, mime=True)

with open(path_fasta, "rb") as f:
    gz = ""
    if path_fasta.endswith(".gz"):
        gz = ".gz"
        f = gzip.GzipFile(fileobj=f, mode='rb')

    f = io.TextIOWrapper(f)

    buffer = f.read(4096)
    seq_type = check_format.from_buffer(buffer).lower()
    f.detach()

    if seq_type != format_to_check:
        raise ValueError(f"Input file ({path_fasta}) does not look like a {format_to_check}")
