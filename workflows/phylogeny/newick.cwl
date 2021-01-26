#!/usr/bin/env cwl-runner

cwlVersion: v1.1

class: CommandLineTool
baseCommand: python

inputs:
    main_py_script:
        type: File
        inputBinding:
            position: 1

    newick:
        type: string
        inputBinding:
            position: 2

    newick_dimensions:
        type: int
        inputBinding:
            position: 3
            prefix: -d

    newick_coverage_matrix:
        type: File
        inputBinding:
            position: 3
    
    newick_metadata:
        type: string
        inputBinding:
            position: 4

    newick_tree:
        type: string
        inputBinding:
            position: 5

outputs:
    metadata_out:
        type: File
        outputBinding:
            glob: 'metadata.tsv'

    newick_tree_out:
        type: File
        outputBinding:
            glob: '*.nwk'