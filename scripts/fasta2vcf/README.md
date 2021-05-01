# Example
```
bash fasta2vcf.sh resources/NC_045512.2.fasta MZ026486.1.fasta MZ026486.1 resources/ensembl-export.csv
```

```
zcat MZ026486.1.vcf.gz

##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=NC_045512.2,length=29903,assembly=NC_045512.2>
##INFO=<ID=ANN,Number=.,Type=String,Description="Associated phenotype">
##bcftools_normVersion=1.10.2+htslib-1.10.2
##bcftools_normCommand=norm -f resources/NC_045512.2.fasta -Ou MZ026486.1.vcf; Date=Sat May  1 12:31:08 2021
##bcftools_annotateVersion=1.10.2+htslib-1.10.2
##bcftools_annotateCommand=annotate --set-id %CHROM\_%POS\_%REF\_%FIRST_ALT -Ov -o -; Date=Sat May  1 12:31:08 2021
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  MZ026486.1
NC_045512.2     174     NC_045512.2_174_G_T     G       T       .       PASS    .       GT      1
NC_045512.2     241     NC_045512.2_241_C_T     C       T       .       PASS    .       GT      1
NC_045512.2     1059    NC_045512.2_1059_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     3037    NC_045512.2_3037_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     4897    NC_045512.2_4897_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     5230    NC_045512.2_5230_G_T    G       T       .       PASS    .       GT      1
NC_045512.2     5512    NC_045512.2_5512_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     6310    NC_045512.2_6310_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     6471    NC_045512.2_6471_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     8326    NC_045512.2_8326_C_T    C       T       .       PASS    .       GT      1
NC_045512.2     10323   NC_045512.2_10323_A_G   A       G       .       PASS    .       GT      1
NC_045512.2     11287   NC_045512.2_11287_GTCTGGTTTT_G  GTCTGGTTTT      G       .       PASS    .       GT      1
NC_045512.2     14408   NC_045512.2_14408_C_T   C       T       .       PASS    .       GT      1
NC_045512.2     16474   NC_045512.2_16474_A_G   A       G       .       PASS    .       GT      1
NC_045512.2     17999   NC_045512.2_17999_C_T   C       T       .       PASS    .       GT      1
NC_045512.2     18657   NC_045512.2_18657_C_T   C       T       .       PASS    .       GT      1
NC_045512.2     21801   NC_045512.2_21801_A_C   A       C       .       PASS    .       GT      1
NC_045512.2     22206   NC_045512.2_22206_A_G   A       G       .       PASS    .       GT      1
NC_045512.2     22280   NC_045512.2_22280_ACTTTACTTG_A  ACTTTACTTG      A       .       PASS    .       GT      1
NC_045512.2     22813   NC_045512.2_22813_G_T   G       T       .       PASS    .       GT      1
NC_045512.2     23012   NC_045512.2_23012_G_A   G       A       .       PASS    .       GT      1
NC_045512.2     23063   NC_045512.2_23063_A_T   A       T       .       PASS    ANN=Fast growing lineage,Increased binding affinity to hACE2 receptor   GT      1
NC_045512.2     23403   NC_045512.2_23403_A_G   A       G       .       PASS    ANN=Moderate effect on transmissibility GT      1
NC_045512.2     23664   NC_045512.2_23664_C_T   C       T       .       PASS    .       GT      1
NC_045512.2     24992   NC_045512.2_24992_G_C   G       C       .       PASS    .       GT      1
NC_045512.2     25563   NC_045512.2_25563_G_T   G       T       .       PASS    .       GT      1
NC_045512.2     25904   NC_045512.2_25904_C_T   C       T       .       PASS    .       GT      1
NC_045512.2     26456   NC_045512.2_26456_C_T   C       T       .       PASS    .       GT      1
NC_045512.2     27676   NC_045512.2_27676_GAACTTTACTCTCCAA_G    GAACTTTACTCTCCAA        G       .       PASS    .       GT      1
NC_045512.2     28253   NC_045512.2_28253_C_T   C       T       .       PASS    .       GT      1
NC_045512.2     28254   NC_045512.2_28254_A_C   A       C       .       PASS    .       GT      1
NC_045512.2     28887   NC_045512.2_28887_C_T   C       T       .       PASS    .       GT      1
NC_045512.2     29856   NC_045512.2_29856_TCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA_T    TCTTAGGAGAATGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA        T       .       PASS    .       GT      1
```


Annotations were downloaded from [here](https://covid-19.ensembl.org/Sars_cov_2/Phenotype/Locations?oa=MONDO:0100096).
