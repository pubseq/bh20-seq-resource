cwlVersion: v1.1
class: CommandLineTool
inputs:
  fasta: File
  metadata: File
outputs: []
requirements:
  DockerRequirement:
    dockerPull: bh20-seq-uploader/import
  NetworkAccess:
    networkAccess: true
baseCommand: bh20-seq-uploader
arguments: [--skip-qc, $(inputs.metadata), $(inputs.fasta)]
