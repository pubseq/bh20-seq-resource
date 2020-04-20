#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
baseCommand: bcftools
arguments:
  - view
  - --no-version
  - -Ou
  - -e'type=ref'
  - --threads=$(inputs.threads)
  - $(inputs.vcf)
inputs:
  - id: vcf
    type: File
  - id: threads
    type: int
outputs:
  - id: bcf
    type: stdout
stdout: $(inputs.vcf.nameroot).without-ref.bcf
