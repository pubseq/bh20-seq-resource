#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/bcftools:1.10.2--hd2cd319_0"
baseCommand: bcftools
arguments:
  - norm
  - -Ob
  - -f
  - $(inputs.ref_fasta)
  - -o
  - $(inputs.output_name)
  - --threads
  - $(inputs.threads)
  - $(inputs.bcf)
inputs:
  - id: ref_fasta
    type: File
  - id: output_name
    type: string
    default: "normalized.bcf"
  - id: threads
    type: int
  - id: bcf
    type: File
outputs:
 - id: normalized_bcf
   type: File
   outputBinding:
     glob: "$(inputs.output_name)"
