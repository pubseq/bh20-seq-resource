#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
baseCommand: bcftools
arguments:
  - view
  - -i
  - 'QUAL>1 && (GT="AA" || GT="Aa")'
  - -Oz
  - --threads=$(inputs.threads)
  - $(inputs.bcf)
inputs:
  - id: threads
    type: int
  - id: bcf
    type: File
    secondaryFiles: [.csi]
outputs:
  - id: vcf
    type: stdout
stdout: out.changes.vcf.gz
