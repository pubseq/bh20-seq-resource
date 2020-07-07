cwlVersion: v1.1
class: CommandLineTool
inputs:
  importScript:
    type: File
    default:
      class: File
      location: download_genbank_data/from_genbank_to_fasta_and_yaml.py
    inputBinding: {position: 1}
  dict:
    type: Directory
    inputBinding:
      prefix: --dict-ontology
      position: 2
    default:
      class: Directory
      location: dict_ontology_standardization
  existing_metadata_from_nuccore:
    type: Directory?
    inputBinding:
      valueFrom: "--skip-request"
      position: 3
outputs:
  fasta_and_yaml:
    type: Directory
    outputBinding:
      glob: fasta_and_yaml
  metadata_from_nuccore:
    type: Directory
    outputBinding:
      glob: metadata_from_nuccore
  accessions:
    type: File?
    outputBinding:
      glob: "*.acc"
  missing_terms:
    type: File
    outputBinding:
      glob: missing_terms.tsv
requirements:
  InitialWorkDirRequirement:
    listing:
      - entry: $(inputs.existing_metadata_from_nuccore)
        entryname: metadata_from_nuccore
  DockerRequirement:
    dockerPull: bh20-seq-uploader/import
  NetworkAccess:
    networkAccess: true
baseCommand: python3
