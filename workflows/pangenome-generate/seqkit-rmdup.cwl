cwlVersion: v1.1
class: CommandLineTool
inputs:
  readsFA: File[]
outputs:
  readsMergeDedup:
    type: File
    outputBinding:
      glob: readsMergeDedup.fasta
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: "quay.io/biocontainers/seqkit:0.7.1--0"
  ResourceRequirement:
    coresMin: 8
    coresMax: 32
    ramMin: $(7 * 1024)
    outdirMin: |
      ${
        var sum = 0;
        for (var i = 0; i < inputs.readsFA.length; i++) {
          sum += inputs.readsFA[i].size;
        }
        return (sum/(1024*1024*1024)+1) + 20;
      }
baseCommand: seqkit
arguments: [rmdup,
            --by-seq,
            --ignore-case,
            -o, readsMergeDedup.fasta,
            $(inputs.readsFA)]
