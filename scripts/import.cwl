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
      location: download_genbank_data/from_genbank_to_fasta_and_yaml.py
    inputBinding: {position: 2}
  dict:
    type: Directory
    default:
      class: Directory
      location: dict_ontology_standardization
    inputBinding: {position: 3}
outputs: []
requirements:
  DockerRequirement:
    dockerPull: bh20-seq-uploader/import
  NetworkAccess:
    networkAccess: true
  WorkReuse:
    enableReuse: false
