import pkg_resources
import tempfile
import magic
import subprocess
import tempfile
import logging
import re
import io
import gzip

log = logging.getLogger(__name__ )

def read_fasta(sequence):
    entries = 0
    bases = []
    label = None
    for line in sequence:
        if line.startswith(">"):
            label = line
            entries += 1
        else:
            bases.append(line)
        if entries > 1:
            log.debug("FASTA file contains multiple entries")
            raise ValueError("FASTA file contains multiple entries")
    return label, bases

def qc_fasta(arg_sequence):
    log.debug("Starting qc_fasta")
    schema_resource = pkg_resources.resource_stream(__name__, "validation/formats")
    with tempfile.NamedTemporaryFile() as tmp:
        tmp.write(schema_resource.read())
        tmp.flush()
        val = magic.Magic(magic_file=tmp.name,
                          uncompress=False, mime=True)

    gz = ""
    if arg_sequence.name.endswith(".gz"):
        sequence = gzip.GzipFile(fileobj=arg_sequence, mode='rb')
        gz = ".gz"
    else:
        sequence = arg_sequence

    sequence = io.TextIOWrapper(sequence)
    r = sequence.read(4096)
    sequence.seek(0)

    seqlabel = r[1:r.index("\n")]
    seq_type = val.from_buffer(r).lower()

    if seq_type == "text/fasta":
        # ensure that contains only one entry
        submitlabel, submitseq = read_fasta(sequence)
        sequence.seek(0)
        sequence.detach()

        with tempfile.NamedTemporaryFile() as tmp1:
            refstring = pkg_resources.resource_string(__name__, "SARS-CoV-2-reference.fasta")
            tmp1.write(refstring)
            tmp1.write(submitlabel.encode("utf8"))
            tmp1.write(("".join(submitseq)).encode("utf8"))
            tmp1.flush()
            subbp = 0
            refbp = 0
            similarity = 0
            try:
                cmd = ["clustalw", "-infile="+tmp1.name,
                       "-quicktree", "-iteration=none", "-type=DNA"]
                print("QC checking similarity to reference")
                print(" ".join(cmd))
                result = subprocess.run(cmd, stdout=subprocess.PIPE)
                res = result.stdout.decode("utf-8")
                g1 = re.search(r"^Sequence 1: [^ ]+ +(\d+) bp$", res, flags=re.MULTILINE)
                refbp = float(g1.group(1))
                g2 = re.search(r"^Sequence 2: [^ ]+ +(\d+) bp$", res, flags=re.MULTILINE)
                subbp = float(g2.group(1))
                g3 = re.search(r"^Sequences \(1:2\) Aligned\. Score: (\d+(\.\d+)?)$", res, flags=re.MULTILINE)
                similarity = float(g3.group(1))

                print(g1.group(0))
                print(g2.group(0))
                print(g3.group(0))
            except Exception as e:
                logging.warn("Error trying to QC against reference sequence using 'clustalw': %s", e)

            if refbp and (subbp/refbp) < .7:
                raise ValueError("QC fail: submit sequence length is shorter than 70% reference")
            if refbp and (subbp/refbp) > 1.3:
                raise ValueError("QC fail: submit sequence length is greater than 130% reference")
            if similarity and similarity < 70.0:
                raise ValueError("QC fail: submit similarity is less than 70%")
            if refbp == 0 or similarity == 0:
                raise ValueError("QC fail")

        return ("sequence.fasta"+gz, seqlabel)
    elif seq_type == "text/fastq":
        sequence.seek(0)
        sequence.detach()
        return ("reads.fastq"+gz, seqlabel)
    else:
        raise ValueError("Sequence file does not look like a DNA FASTA or FASTQ")
