#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: CommandLineTool

baseCommand: awk

inputs:
  consensus_regex:
    type: string
    inputBinding:
        position: 1

  coverage_tsv:
    type: File
    inputBinding:
        position: 2
 
outputs:
  awk_coverage_matrix:
    type: stdout

stdout: coverage.no_consensus.tsv