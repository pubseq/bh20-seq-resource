import pkg_resources
import tempfile
import magic

def qc_fasta(sequence):
    schema_resource = pkg_resources.resource_stream(__name__, "validation/formats")
    with tempfile.NamedTemporaryFile() as tmp:
        tmp.write(schema_resource.read())
        tmp.flush()
        val = magic.Magic(magic_file=tmp.name,
                          uncompress=False, mime=True)
    seq_type = val.from_buffer(sequence.read(4096)).lower()
    sequence.seek(0)
    if seq_type == "text/fasta":
        # ensure that contains only one entry
        entries = 0
        for line in sequence:
            if line.startswith(">"):
                entries += 1
            if entries > 1:
                raise ValueError("FASTA file contains multiple entries")
                break
        sequence.seek(0)
        return "reads.fastq"
    elif seq_type == "text/fastq":
        return "sequence.fasta"
    else:
        raise ValueError("Sequence file does not look like FASTA or FASTQ")
