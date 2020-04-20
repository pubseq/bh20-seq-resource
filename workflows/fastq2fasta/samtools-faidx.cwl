#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
doc: "samtools sort, sort given bam file"
requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input_fasta)
baseCommand: [samtools, faidx]
inputs:
  input_fasta:
    type: File
    label: "Input fasta"
    inputBinding:
      position: 1
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
outputs:
  indexed_fasta:
    type: File
    outputBinding:
      glob: "$(inputs.input_fasta.basename)"
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
      - .fai
  stdout: stdout
  stderr: stderr
stdout: samtools-sort-stdout.log
stderr: samtools-sort-stderr.log
