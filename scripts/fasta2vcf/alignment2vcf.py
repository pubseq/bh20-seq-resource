# From https://github.com/boxiangliu/covseq

from collections import defaultdict

def parse_align(align_fn):
	header_ref = ''
	ref = ''

	header_qry = ''
	qry = ''
	n_record = 0
	with open(align_fn, "r") as f:
		for line in f:
			if line.startswith(">"):
				if n_record == 0:
					is_ref = True
					n_record += 1

					header_ref = line.strip()
				else:
					is_ref = False

					header_qry = line.strip()
			else:
				if is_ref:
					ref += line.strip()
				else:
					qry += line.strip()

	return (header_ref, ref), (header_qry, qry)

def align2variant(ref, qry):
	assert len(ref) == len(qry)
	ref_coord = 0
	qry_coord = 0
	ref_variant = defaultdict(str)
	qry_variant = defaultdict(str)
	r0 = ""
	q0 = ""

	for r, q in zip(ref,qry):
		r = r.upper().replace("U", "T")
		q = q.upper().replace("U", "T")

		if r != "-":
			ref_coord += 1
		if q != "-":
			qry_coord = ref_coord
		if r != "-" and q != "-":
			r0 = r
			q0 = q

		if r == q:
			pass
		elif r == "n" or q == "n":
			pass
		elif r == "-":
			if ref_variant[ref_coord] == "":
				ref_variant[ref_coord] = r0
				qry_variant[ref_coord] = q0
			ref_variant[ref_coord] += r
			qry_variant[ref_coord] += q
		elif q == "-":
			if ref_variant[qry_coord] == "":
				ref_variant[qry_coord] = r0
				qry_variant[qry_coord] = q0
			ref_variant[qry_coord] += r
			qry_variant[qry_coord] += q
		elif r != q:
			ref_variant[ref_coord] = r
			qry_variant[ref_coord] = q
		else:
			raise Exception("Error!")
	return ref_variant, qry_variant

def save_vcf(ref_variant, qry_variant, qry_name, out_fn):
	assert len(ref_variant) == len(qry_variant)

	with open(out_fn, "w") as f:
		f.write('##fileformat=VCFv4.2\n')
		f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
		f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
		f.write('##contig=<ID=NC_045512.2,length=29903,assembly=NC_045512.2>\n')
		f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{qry_name}\n")
		for coord in ref_variant.keys():
			if coord == 0: # skip coord = 0
				continue

			rv = ref_variant[coord].replace("-","")
			qv = qry_variant[coord].replace("-","")

			rv_is_canonical = all([x in ["A","T","G","C"] for x in rv])
			qv_is_canonical = all([x in ["A","T","G","C"] for x in qv])

			if rv_is_canonical and qv_is_canonical:
				f.write(f"NC_045512.2\t{coord}\t.\t{rv}\t{qv}\t.\tPASS\t.\tGT\t1\n")


def filter_polya(vcf_fn):
	with open(vcf_fn, "r") as fin, open(f"{vcf_fn}.tmp", "w") as fout:
		for line in fin:
			if line.startswith("#"):
				fout.write(line)
			else:
				split_line = line.strip().split("\t")
				pos = split_line[1]
				ref = split_line[3]
				if pos == "29870" and ref.endswith("AAAAAAAAAA"):
					pass
				else:
					fout.write(line)
	os.remove(vcf_fn)
	os.rename(f"{vcf_fn}.tmp", vcf_fn)


import sys
import os

path_reference = sys.argv[1]
path_alignment = sys.argv[2]
output_prefix = sys.argv[3]

print(path_reference)
print(path_alignment)
print("Parsing alignment to get reference and query sequences")
(header_ref, ref), (header_qry, qry) = parse_align(path_alignment)

print("Identifying variants")
ref_variant, qry_variant = align2variant(ref, qry)

print(f"Saving variants to VCF")
out_fn = f'{output_prefix}.vcf'
save_vcf(ref_variant, qry_variant, output_prefix, out_fn)

print("Removing variants in the Poly-A tail.")
filter_polya(out_fn)

#print(f"Normalize, update ID, and index.")
#postprocess_vcf(out_fn, path_reference)
