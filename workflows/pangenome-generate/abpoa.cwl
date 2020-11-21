cwlVersion: v1.1
class: CommandLineTool
inputs:
  readsFA: File
  script:
    type: File
    default: {class: File, location: relabel-seqs.py}
outputs:
  abpoaGFA:
    type: stdout
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/abpoa:1.0.5--hed695b0_0"
  ResourceRequirement:
    coresMin: 1
    ramMin: $(15 * 1024)
    outdirMin: $(Math.ceil(inputs.readsFA.size/(1024*1024*1024) + 20))
baseCommand: abpoa
stdout: $(inputs.readsFA.nameroot).O0.gfa
arguments: [
    $(inputs.readsFA),
    -r 3,
    -O, '0'
]
