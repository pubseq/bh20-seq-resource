#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: CommandLineTool
doc: "Workflow to go from YAML (metadata) + FASTA (sequence) to TTL (metadata)"

inputs:
  path_fasta:
    type: File
    inputBinding:
      position: 1
  path_yaml:
    type: File
    inputBinding:
      position: 2

steps:
  check_format:
    in: {path_fasta: path_fasta, path_valid_formats: '../../bh20sequploader/validation/formats', format_to_check: 'text/fasta'}
    #out: true/false or nothing and it has to block the execution if the format is wrong
    run: check_format.cwl
  check_metadata:
    # input and output
    # run: check_metadata.cwl
  check_header:
    # id_fasta has to be equal to id_yaml
    # run: check_header.cwl
  check_sequence:
    # The sequence has to be similar to the reference
    # run: check_sequence.cwl

