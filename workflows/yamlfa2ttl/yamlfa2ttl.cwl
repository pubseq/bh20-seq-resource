#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: Workflow
doc: "Workflow to go from YAML (metadata) + FASTA (sequence) to TTL (metadata)"

inputs:
  path_fasta:
    type: string
    doc: input fasta to validate

  format_to_check:
    type: string
    default: text/fasta

steps:
  check_format:
    in:
      path_fasta: path_fasta
      format_to_check: format_to_check
    out: []
    run: check_format.cwl

outputs: []
