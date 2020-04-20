#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
baseCommand: bcftools
arguments:
  - concat
  - -Ou
  - -o
  - $(inputs.output_name)
  - $(inputs.bcf_files)
inputs:
  - id: output_name
    type: string
    default: "merged.bcf"
  - id: bcf_files
    type: File[]
outputs:
  - id: merged_bcf
    type: File
    outputBinding:
      glob: "$(inputs.output_name)"
