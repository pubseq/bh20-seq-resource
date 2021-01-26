#!/usr/bin/env cwl-runner

cwlVersion: v1.1

class: CommandLineTool
baseCommand: odgi

inputs:
  odgi_paths:
    type: string
    inputBinding:
      position: 1

  odgi_graph:
    type: File
    inputBinding:
        position: 2
        prefix: -i

  haplotypes:
    type: boolean
    inputBinding:
      position: 4
      prefix: -H

  threads:
    type: int
    inputBinding:
        position: 5
        prefix: -t

outputs:
  coverage_matrix:
    type: stdout

stdout: coverage.tsv