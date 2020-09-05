#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
  ShellCommandRequirement: {}
baseCommand: bcftools
arguments:
  - consensus
  - -i
  - 'QUAL > 10 && GT="a"'
  - -Hla
  - -f
  - $(inputs.ref_fasta)
  - $(inputs.vcf)
  - {shellQuote: false, valueFrom: "|"}
  - sed
  - "s/^>.*/>$(inputs.sample_id)/g"
inputs:
  - id: ref_fasta
    type: File
  - id: vcf
    type: File
    secondaryFiles: [.csi]
  - id: sample_id
    type: string
outputs:
  - id: out_fasta
    type: stdout
stdout: sequence.fasta
