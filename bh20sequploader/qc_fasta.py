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

def qc_fasta(arg_sequence, check_with_mimimap2=True):
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

        if check_with_mimimap2:
            with tempfile.NamedTemporaryFile() as tmp1:
                with tempfile.NamedTemporaryFile() as tmp2:
                    refstring = pkg_resources.resource_string(__name__, "SARS-CoV-2-reference.fasta")
                    tmp1.write(refstring)
                    tmp1.flush()
                    tmp2.write(submitlabel.encode("utf8"))
                    tmp2.write(("".join(submitseq)).encode("utf8"))
                    tmp2.flush()

                    similarity = 0
                    try:
                        cmd = ["minimap2", "-c -x asm20", tmp1.name, tmp2.name]
                        logging.info("QC checking similarity to reference")
                        logging.info(" ".join(cmd))
                        result = subprocess.run(cmd, stdout=subprocess.PIPE)
                        result.check_returncode()
                        res = result.stdout.decode("utf-8")
                        mm = res.split("\t")
                        if len(mm) >= 10:
                            # divide Number of matching bases in the mapping / Target sequence length
                            similarity = (float(mm[9]) / float(mm[6])) * 100.0
                        else:
                            similarity = 0
                    except Exception as e:
                        logging.warn("QC against reference sequence using 'minimap2': %s", e, exc_info=e)

                    if similarity < 70.0:
                        raise ValueError(
                            "QC fail for {}: alignment to reference was less than 70%% (was %2.2f%%)".format(
                                seqlabel, similarity
                            ))

        return "sequence.fasta" + gz, seqlabel, seq_type
    elif seq_type == "text/fastq":
        sequence.seek(0)
        sequence.detach()
        return "reads.fastq" + gz, seqlabel, seq_type
    else:
        raise ValueError("Sequence file ({}) does not look like a DNA FASTA or FASTQ".format(arg_sequence))
