#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
doc: "samtools sort, sort given bam file"
requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/samtools:1.9--h8571acd_11
baseCommand: [samtools, sort]
inputs:
  threads:
    type: int
    default: 4
    inputBinding:
      prefix: -@
  tmpfile:
    type: string
    default: sort.tmp
    label: "Write temporary files to PREFIX.nnnn.bam"
    inputBinding:
      prefix: -T
  output_bam:
    type: string
    default: aln.sorted.bam
    label: "Write final output to FILENAME"
    inputBinding:
      prefix: -o
  input_bamfile:
    type: File
    label: "Input bamfile"
    inputBinding:
      position: 1

outputs:
  sorted_bam:
    type: File
    outputBinding:
      glob: "$(inputs.output_bam)"
  stdout: stdout
  stderr: stderr
stdout: samtools-sort-stdout.log
stderr: samtools-sort-stderr.log
