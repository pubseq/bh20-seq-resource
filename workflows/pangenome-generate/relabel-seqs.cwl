cwlVersion: v1.1
class: CommandLineTool
inputs:
  readsFA: File[]
  subjects: string[]
outputs:
  relabeledSeqs:
    type: stdout
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entry: {$include: relabel-seqs.py}
        entryname: relabel-seqs.py
hints:
  DockerRequirement:
    dockerPull: commonworkflowlanguage/cwltool_module
stdout: relabeledSeqs.fasta
baseCommand: [python, relabel-seqs.py]
