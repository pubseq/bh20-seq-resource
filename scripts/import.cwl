cwlVersion: v1.1
class: CommandLineTool
baseCommand: python3
inputs:
  scripts:
    type: File
    default:
      class: File
      location: import_to_arvados.py
    inputBinding: {position: 1}
  importScript:
    type: File
    default:
      class: File
      location: from_genbank_to_fasta_and_yaml.py
    inputBinding: {position: 2}
outputs: []
requirements:
  DockerRequirement:
    dockerPull: bh20-seq-uploader/import
  NetworkAccess:
    networkAccess: true
  WorkReuse:
    workReuse: false
