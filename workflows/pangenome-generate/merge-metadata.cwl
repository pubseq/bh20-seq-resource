cwlVersion: v1.1
class: CommandLineTool
hints:
  DockerRequirement:
    dockerPull: commonworkflowlanguage/cwltool_module
inputs:
  metadata: File[]
  subjects: string[]
  metadataSchema:
    type: File
    inputBinding: {position: 2}
  originalLabels:
    type: File
    inputBinding: {position: 3}
  dups:
    type: File?
    inputBinding: {position: 4}
  script:
    type: File
    inputBinding: {position: 1}
    default: {class: File, location: merge-metadata.py}
outputs:
  merged: stdout
stdout: mergedmetadata.ttl
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing: |
          ${
          var i = 0;
          var b = 1;
          var out = [];
          for (; i < inputs.metadata.length; i++) {
            var block = [];
            var sub = [];
            for (; i < (b*150) && i < inputs.metadata.length; i++) {
              block.push(inputs.metadata[i]);
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
baseCommand: python
