cwlVersion: v1.1
class: CommandLineTool
inputs:
  inputGFA: File
outputs:
  odgiGraph:
    type: File
    outputBinding:
      glob: $(inputs.inputGFA.nameroot).unchop.sorted.odgi
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "odgi-bash-binutils:latest"
  ResourceRequirement:
    coresMin: 4
    ramMin: $(15 * 1024)
    outdirMin: $(Math.ceil((inputs.inputGFA.size/(1024*1024*1024)+1) * 2))
  InitialWorkDirRequirement:
    # Will fail if input file is not writable (odgi bug)
    listing:
      - entry: $(inputs.inputGFA)
        writable: true
arguments:
  - "sh"
  - "-c"
  - >-
    odgi build -g '$(inputs.inputGFA.path)' -o - | odgi unchop -i - -o - |
    odgi sort -i - -p s -o $(inputs.inputGFA.nameroot).unchop.sorted.odgi
