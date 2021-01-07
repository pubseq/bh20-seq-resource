~/.config/guix/current/bin/guix environment -C guix --ad-hoc cwltool python python-biopython python-requests python-dateutil python-magic ruby
cwltool --preserve-environment PYTHONPATH yamlfa2ttl.cwl --path_fasta ~/bh20-seq-resource/example/sequence.fasta

cwltool --no-container --preserve-environment GUIX_ENVIRONMENT --preserve-environment PYTHONPATH yamlfa2ttl.cwl --path_fasta ~/bh20-seq-resource/example/sequence.fasta



#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: Workflow
doc: "Workflow to go from YAML (metadata) + FASTA (sequence) to TTL (metadata)"

inputs:
  path_fasta:
    type: string
    doc: input FASTA to validate

  format_to_check:
    type: string
    default: text/fasta

  path_yaml:
    type: string
    doc: input YAML to validate and convert in TTL

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

  check_metadata:
    in:
      path_yaml: path_yaml
    doc: the input metadata information to put in the knowledge graph
    out: []
    run: check_metadata.cwl

outputs: []
