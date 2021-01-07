#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: CommandLineTool
baseCommand: python3

inputs:
  script:
    type: File
    inputBinding: {position: 1}
    default: {class: File, location: check_sequence.py}
  path_fasta:
    type: string
    inputBinding: {position: 2}
  path_sars_cov_2_reference_fasta:
    type: File
    inputBinding: {position: 3}
    default: {class: File, location: ../../bh20sequploader/SARS-CoV-2-reference.fasta}

outputs: []
