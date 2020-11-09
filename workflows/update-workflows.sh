#!/bin/sh
arvados-cwl-runner --project-uuid=lugli-j7d0g-5hswinmpyho8dju --update-workflow=lugli-7fd4e-2zp9q4jo5xpif9y fastq2fasta/fastq2fasta.cwl
arvados-cwl-runner --project-uuid=lugli-j7d0g-5hswinmpyho8dju --update-workflow=lugli-7fd4e-mqfu9y3ofnpnho1 pangenome-generate/collect-seqs.cwl
