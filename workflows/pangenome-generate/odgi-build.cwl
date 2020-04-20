cwlVersion: v1.1
class: CommandLineTool
inputs:
  inputGFA: File
outputs:
  odgiGraph:
    type: File
    outputBinding:
      glob: $(inputs.inputGFA.nameroot).odgi
requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/odgi:v0.3--py37h8b12597_0"
  ResourceRequirement:
    coresMin: 4
    ramMin: $(7 * 1024)
    outdirMin: $(Math.ceil((inputs.inputGFA.size/(1024*1024*1024)+1) * 2))
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.inputGFA)
        writable: true
arguments: [odgi, build, -g, $(inputs.inputGFA), -s, -o, -,
            {shellQuote: false, valueFrom: "|"},
            odgi, sort, -i, -, -p, s, -o, $(inputs.inputGFA.nameroot).odgi]
