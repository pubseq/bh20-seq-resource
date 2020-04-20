#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
doc: "samtools view to convert sam format to bam format"
requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11
baseCommand: [samtools, view]
inputs:
  threads:
    type: int
    label: "Number of additional threads to use"
    default: 4
    inputBinding:
      prefix: -@
  output_bam:
    type: boolean
    label: "output BAM"
    default: true
    inputBinding:
      prefix: -b
  output_filename:
    type: string
    label: "output file name"
    default: "aln.bam"
    inputBinding:
      prefix: -o
  input_file:
    type: File
    label: "input file"
    inputBinding:
      position: 1
  include_header:
    type: boolean
    label: "include the header in the output"
    default: false
    inputBinding:
      prefix: -h
  ignore_previous_version:
    type: boolean
    label: "ignored for compatibility with previous samtools versions"
    default: false
    inputBinding:
      prefix: -S
  filter_alignments:
    type: string?
    label: "Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0]."
    inputBinding:
      prefix: -F
  skip_alignments:
    type: int?
    label: "Skip alignments with MAPQ smaller than INT [0]."
    inputBinding:
      prefix: -q
outputs:
  bam:
    type: File
    outputBinding:
      glob: "$(inputs.output_filename)"
  stdout: stdout
  stderr: stderr
stdout: samtools-view-stdout.log
stderr: samtools-view-stderr.log
