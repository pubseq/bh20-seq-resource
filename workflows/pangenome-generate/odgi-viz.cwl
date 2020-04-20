cwlVersion: v1.1
class: CommandLineTool
inputs:
  inputODGI: File
outputs:
  odgiPNG:
    type: File
    outputBinding:
      glob: $(inputs.inputODGI.nameroot).png
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/odgi:v0.3--py37h8b12597_0"
  ResourceRequirement:
    coresMin: 4
    ramMin: $(7 * 1024)
    outdirMin: 1
baseCommand: [odgi, viz]
arguments: [-i, $(inputs.inputODGI),
            -o, $(inputs.inputODGI.nameroot).png,
            -x, "50000",
            -y, "500",
            -R,
            -P, "4"]
