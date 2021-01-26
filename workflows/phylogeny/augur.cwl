#!/usr/bin/env cwl-runner

cwlVersion: v1.1

class: CommandLineTool
baseCommand: bash

requirements:
  InitialWorkDirRequirement:
    listing:
      - $(inputs.dataDir)

inputs:
    nextstrain_bash_script:
        type: File
        inputBinding:
            position: 1

    newick_tree_2:
        type: File
        inputBinding:
            position: 2

    metadata_newick:
        type: File
        inputBinding:
            position: 3

    dataDir:
        type: Directory

outputs:
    newick_json:
        type: File
        outputBinding:
            glob: 'covid.json'