#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: CommandLineTool
baseCommand: python3

inputs:
  script:
    type: File
    inputBinding: {position: 1}
    default: {class: File, location: check_metadata.py}
  path_yaml:
    type: string
    inputBinding: {position: 2}
  path_schema_yaml:
    type: File
    inputBinding: {position: 3}
    default: {class: File, location: ../../bh20sequploader/bh20seq-schema.yml}
  path_shex_rdf:
    type: File
    inputBinding: {position: 4}
    default: {class: File, location: ../../bh20sequploader/bh20seq-shex.rdf}

outputs: []
