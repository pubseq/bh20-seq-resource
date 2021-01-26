#!/usr/bin/env cwl-runner

cwlVersion: v1.1

class: CommandLineTool
baseCommand: python

inputs:
    main_py_script:
        type: File
        inputBinding:
            position: 1

    metadata:
        type: string
        inputBinding:
            position: 2

    coverage_matrix:
        type: File
        inputBinding:
            position: 3

    coverage_matrix_with_metadata:
        type: string
        inputBinding:
            position: 4

outputs:
    coverage_matrix_with_metadata_out:
        type: File
        outputBinding:
            glob: '*.metadata.tsv'