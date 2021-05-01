# Example
```
bash fasta2vcf.sh resources/NC_045512.2.fasta MZ026486.1.fasta MZ026486.1
```

```
zcat MZ026486.1.vcf.gz | head -n 15

##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=NC_045512.2,length=29903,assembly=NC_045512.2>
##bcftools_normVersion=1.10.2+htslib-1.10.2
##bcftools_normCommand=norm -f resources/NC_045512.2.fasta -Ou MZ026486.1.vcf; Date=Sat May  1 11:01:01 2021
##bcftools_annotateVersion=1.10.2+htslib-1.10.2
##bcftools_annotateCommand=annotate --set-id %CHROM\_%POS\_%REF\_%FIRST_ALT -Ov -o -; Date=Sat May  1 11:01:01 2021
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  MZ026486.1
NC_045512.2     174     NC_045512.2_174_G_T     G       T       .       PASS    .       GT      1
NC_045512.2     241     NC_045512.2_241_C_T     C       T       .       PASS    .       GT      1
NC_045512.2     1059    NC_045512.2_1059_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     3037    NC_045512.2_3037_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     4897    NC_045512.2_4897_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     5230    NC_045512.2_5230_G_T    G       T       .       PASS    .       GT      1
```