import sys
import os

path_vcf = sys.argv[1]
path_annotation = sys.argv[2]

variant_to_phenotypes_dict = {}

with open(path_annotation) as f:
    f.readline()  # Skip header

    for line in f:
        if line.startswith('MN908947'):
            variant = line.split(',')[0]
            pos, ref, alt = variant.split(':')[1:]

            # Ugly, Pjotr will not like it
            f.readline()
            f.readline()
            phenotype = f.readline().split(',')[-3].strip('"')

            if (pos, ref, alt) not in variant_to_phenotypes_dict:
                variant_to_phenotypes_dict[(pos, ref, alt)] = []
            variant_to_phenotypes_dict[(pos, ref, alt)].append(phenotype)

new_row_in_header = '##INFO=<ID=ANN,Number=.,Type=String,Description="Associated phenotype">\n'

with open(path_vcf) as fin, open(f"{path_vcf}.tmp", "w") as fout:
    for line in fin:
        if line:
            if line.startswith('#CHROM'):
                if new_row_in_header:
                    fout.write(new_row_in_header)
                    new_row_in_header = ''

            if not line.startswith("#"):
                split_line = line.strip().split("\t")
                pos = split_line[1]
                ref, alt = split_line[3:5]

                if (pos, ref, alt) in variant_to_phenotypes_dict:
                    split_line[7] = "ANN={}".format(
                        ",".join(variant_to_phenotypes_dict[(pos, ref, alt)])
                    )
                line = "\t".join(split_line) + "\n"

        fout.write(line)
os.remove(path_vcf)
os.rename(f"{path_vcf}.tmp", path_vcf)
