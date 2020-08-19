cwlVersion: v1.1
class: CommandLineTool
baseCommand: python
inputs:
  script:
    type: File
    default:
      class: File
      location: dups2metadata.py
    inputBinding: {position: 1}
  metadata:
    type: File
    inputBinding: {position: 2}
  dups:
    type: File?
    inputBinding: {position: 3}
stdout: mergedmetadata.ttl
outputs:
  merged: stdout
