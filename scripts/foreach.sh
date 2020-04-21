#!/bin/sh
rm -rf validated fasta_and_yaml_*
mkdir -p validated
./from_genbank_to_fasta_and_yaml.py
fasta_files=$(find fasta_and_yaml_20200421/ -name "*.fasta")
for f in $fasta_files ; do
    yaml=$(echo $f | rev | cut -c7- | rev).yaml
    echo $f
    echo $yaml
    if bh20-seq-uploader --validate $f $yaml ; then
	sz=$(stat --format=%s $f)
	if test $sz -gt 20000 ; then
	    mv $f $yaml validated
	else
	    echo "Fasta file too small"
	fi
    fi
done
