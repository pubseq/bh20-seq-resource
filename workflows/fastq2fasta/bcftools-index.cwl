#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
  InitialWorkDirRequirement:
    listing:
      - $(inputs.bcf)
baseCommand: bcftools
arguments:
  - index
  - $(inputs.bcf)
inputs:
  - id: bcf
    type: File
outputs:
  - id: indexed
    type: File
    outputBinding:
      glob: "$(inputs.bcf.basename)"
    secondaryFiles:
      - .csi
