#!/usr/bin/env bash

path_reference=$1
path_query=$2
output_prefix=$3

echo "Contatenating reference and query in the same FASTA file"
cat $path_reference $path_query > ref+qry.fasta

echo "Aligning reference and query with MAFFT"
mafft ref+qry.fasta > ref+qry.alignment

python3 alignment2vcf.py $path_reference ref+qry.alignment $output_prefix

java -jar /home/tools/snpEff/5.0e/snpEff.jar NC_045512.2 $output_prefix.vcf | bgzip -c > $output_prefix.annotated.vcf.gz && tabix -p vcf $output_prefix.annotated.vcf.gz

echo "Removing temporary files"
rm snpEff_genes.txt snpEff_summary.html ref+qry.fasta ref+qry.alignment $output_prefix.vcf