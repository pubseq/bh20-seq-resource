#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: CommandLineTool
doc: string
requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bwa:0.7.17--h84994c4_5
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input_fasta)
baseCommand: [bwa, index]
inputs:
  input_fasta:
    type: File
    label: "input fasta file"
    inputBinding:
      position: 1
outputs:
  indexed_fasta:
    type: File
    outputBinding:
      glob: $(inputs.input_fasta.basename)
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
  stdout: stdout
  stderr: stderr
stdout: bwa-index-stdout.log
stderr: bwa-index-stderr.log
