cwlVersion: v1.1
class: CommandLineTool
inputs:
  readsFA: File
  script:
    type: File
    default: {class: File, location: relabel-seqs.py}
outputs:
  spoaGFA:
    type: stdout
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/spoa:3.0.2--hc9558a2_0"
  ResourceRequirement:
    coresMin: 1
    ramMin: $(15 * 1024)
    outdirMin: $(Math.ceil(inputs.readsFA.size/(1024*1024*1024) + 20))
baseCommand: spoa
stdout: $(inputs.readsFA.nameroot).g6.gfa
arguments: [
    $(inputs.readsFA),
    -G,
    -g, '-6'
]
