cwlVersion: v1.1
class: CommandLineTool
inputs:
  readsFA: File[]
  subjects: string[]
outputs:
  relabeledSeqs:
    type: File
    outputBinding:
      glob: relabeledSeqs.fasta
  originalLabels:
    type: File
    outputBinding:
      glob: originalLabels.ttl
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: {$include: relabel-seqs.py}
        entryname: relabel-seqs.py
hints:
  DockerRequirement:
    dockerPull: commonworkflowlanguage/cwltool_module
stdout:
baseCommand: [python, relabel-seqs.py]
