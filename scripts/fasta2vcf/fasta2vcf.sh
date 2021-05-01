#!/usr/bin/env bash

path_reference=$1
path_query=$2
output_prefix=$3
path_annotation=$4

echo "Contatenating reference and query in the same FASTA file"
cat $path_reference $path_query > ref+qry.fasta

echo "Aligning reference and query with MAFFT"
mafft ref+qry.fasta > ref+qry.alignment

python3 alignment2vcf.py $path_reference ref+qry.alignment $output_prefix

python3 simpleVcfAnnotation.py $output_prefix.vcf $path_annotation

bcftools norm -f $path_reference $output_prefix.vcf -Ou | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -Ov -o - | bgzip -c > $output_prefix.vcf.gz
#tabix -p vcf $output_prefix.vcf.gz

#java -jar /home/tools/snpEff/5.0e/snpEff.jar NC_045512.2 $output_prefix.vcf | bgzip -c > $output_prefix.annotated.vcf.gz && tabix -p vcf $output_prefix.annotated.vcf.gz

echo "Removing temporary files"
#rm snpEff_genes.txt snpEff_summary.html
rm ref+qry.fasta ref+qry.alignment $output_prefix.vcf