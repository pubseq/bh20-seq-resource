#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: CommandLineTool
baseCommand: python3
inputs:
  script:
    type: File
    inputBinding: {position: 1}
    default: {class: File, location: check_format.py}
  path_fasta:
    type: string
    inputBinding: {position: 2}
  format_to_check:
    type: string
    inputBinding: {position: 3}
  path_valid_formats:
    type: File
    inputBinding: {position: 4}
    default: {class: File, location: ../../bh20sequploader/validation/formats}
outputs: []

