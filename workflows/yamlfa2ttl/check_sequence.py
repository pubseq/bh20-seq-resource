import subprocess
import tempfile
import sys

path_fasta = sys.argv[1]
path_sars_cov_2_reference_fasta = sys.argv[2]


def read_single_fasta(path_fasta):
    with open(path_fasta) as f:
        entries = 0
        header = None
        sequence = []
        for line in f:
            line = line.strip()

            if line.startswith(">"):
                header = line
                entries += 1
            else:
                sequence.append(line)

            if entries > 1:
                raise ValueError("Input FASTA contains multiple entries")

    return header, ''.join(sequence)


print("FASTA QC: checking similarity to the reference", file=sys.stderr)

header, sequence = read_single_fasta(path_fasta)

similarity = 0

with tempfile.NamedTemporaryFile() as tmp_fasta:
    with tempfile.NamedTemporaryFile() as tmp_sars_cov_2_reference_fasta:
        with open(path_sars_cov_2_reference_fasta, 'rb') as f:
            tmp_sars_cov_2_reference_fasta.write(f.read())
        tmp_sars_cov_2_reference_fasta.flush()

        tmp_fasta.write(f'>{header}\n'.encode("utf8"))
        tmp_fasta.write(sequence.encode("utf8"))
        tmp_fasta.flush()

        cmd = [
            "minimap2", "-c", "-x", "asm20",
            tmp_sars_cov_2_reference_fasta.name, tmp_fasta.name
        ]
        print(" ".join(cmd), file=sys.stderr)

        result = subprocess.run(cmd, stdout=subprocess.PIPE)
        result.check_returncode()

        paf_output = result.stdout.decode("utf-8")
        paf_output = paf_output.split("\t")
        if len(paf_output) >= 10:
            # Number of matching bases in the mapping / Target sequence length
            similarity = (float(paf_output[9]) / float(paf_output[6])) * 100.0

if similarity < 70.0:
    raise ValueError(
        f"FASTA QC fail for '{header}': the similarity against the reference {similarity} was less than 70%')"
    )
