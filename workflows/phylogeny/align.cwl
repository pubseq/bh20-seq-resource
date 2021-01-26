#!/usr/bin/env cwl-runner

cwlVersion: v1.1

class: CommandLineTool
baseCommand: pggb

inputs:
  threads:
    type: int
    inputBinding:
        position: 1
        prefix: -t

  pggb_wfmash:
    type: boolean
    inputBinding:
        position: 2
        prefix: --wfmash

  pggb_fasta:
    type: File
    inputBinding:
        position: 3
        prefix: -i

  pggb_mash_k_mer:
    type: int
    inputBinding:
        position: 4
        prefix: -K

  pggb_map_percent_identity:
    type: int
    inputBinding:
        position: 5
        prefix: -p

  pggb_num_secondary_mappings:
    type: int
    inputBinding:
        position: 6
        prefix: -n

  pggb_segment_length:
    type: int
    inputBinding:
        position: 7
        prefix: -s

  pggb_output_dir:
    type: string
    inputBinding:
        position: 8
        prefix: -o

outputs:
    pggb_odgi_graph:
        type: File
        outputBinding:
            glob: '*.smooth.og'