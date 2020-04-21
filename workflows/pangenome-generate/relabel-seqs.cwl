cwlVersion: v1.1
class: CommandLineTool
inputs:
  readsFA: File[]
  subjects: string[]
  script:
    type: File
    default: {class: File, location: relabel-seqs.py}
    inputBinding: {}
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
    listing: |
          ${
          var i = 0;
          var b = 1;
          var out = [];
          for (; i < inputs.readsFA.length; i++) {
            var block = [];
            var sub = [];
            for (; i < (b*150) && i < inputs.readsFA.length; i++) {
              block.push(inputs.readsFA[i]);
              sub.push(inputs.subjects[i]);
            }
            out.push({
              entryname: "block"+b,
              entry: JSON.stringify(block)
            });
            out.push({
              entryname: "subs"+b,
              entry: JSON.stringify(sub)
            });
            b++;
          }
          return out;
          }
hints:
  DockerRequirement:
    dockerPull: commonworkflowlanguage/cwltool_module
stdout:
baseCommand: [python]
