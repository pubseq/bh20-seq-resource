cwlVersion: v1.1
class: CommandLineTool
inputs:
  readsFA:
    type: File
    inputBinding: {position: 2}
  script:
    type: File
    inputBinding: {position: 1}
    default: {class: File, location: sort_fasta_by_quality_and_len.py}
stdout: $(inputs.readsFA.nameroot).sorted_by_quality_and_len.fasta
outputs:
  sortedReadsFA:
    type: stdout
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
baseCommand: [python]
