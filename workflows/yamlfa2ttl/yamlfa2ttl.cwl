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
    doc: the input has to be a valid FASTA format file
    out: []
    run: check_format.cwl

  check_sequence:
    in:
      path_fasta: path_fasta
    doc: the input sequence has to be enough similar to the reference
    out: []
    run: check_sequence.cwl

outputs: []
