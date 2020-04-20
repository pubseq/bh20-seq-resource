cwlVersion: v1.1
class: CommandLineTool
inputs:
  readsFA: File
outputs:
  readsPAF: stdout
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/minimap2:2.17--h8b12597_1"
  ResourceRequirement:
    coresMin: 8
    coresMax: 32
    ramMin: $(7 * 1024)
    outdirMin: $(Math.ceil(inputs.readsFA.size/(1024*1024*1024) + 20))
stdout: $(inputs.readsFA.nameroot).paf
baseCommand: minimap2
arguments: [-cx, asm20,
            -w, "1",
            -t, $(runtime.cores),
            $(inputs.readsFA),
            $(inputs.readsFA)]
