#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.1
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/freebayes:1.3.2--py37hc088bd4_0"
baseCommand: freebayes
arguments: [
  --bam, $(inputs.bam),
  # --region=$(inputs.contig):1-$(inputs.contig_end)
  --ploidy, "1",
  -f, $(inputs.ref_fasta)]
inputs:
  - id: bam
    type: File
  # - id: contig
  #   type: string
  # - id: contig_end
  #   type: int
  - id: ref_fasta
    type: File
outputs:
  - id: vcf
    type: stdout
stdout: var.vcf
