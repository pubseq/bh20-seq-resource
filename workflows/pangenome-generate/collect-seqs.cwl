cwlVersion: v1.1
class: CommandLineTool
$namespaces:
  arv: "http://arvados.org/cwl#"
  cwltool: "http://commonwl.org/cwltool#"
requirements:
  arv:APIRequirement: {}
  arv:RuntimeConstraints:
    outputDirType: keep_output_dir
  DockerRequirement:
    dockerPull: arvados/jobs:2.0.3
  WorkReuse:
    enableReuse: false
  ResourceRequirement:
    coresMin: 1
    ramMin: 1024
baseCommand: python3
inputs:
  script:
    type: File
    default:
      class: File
      location: collect-seqs.py
    inputBinding: {position: 1}
  src_project:
    type: string
    inputBinding: {position: 2}
  metadataSchema:
    type: File
    inputBinding: {position: 3}
  exclude:
    type: File?
    inputBinding: {position: 4}
outputs:
  relabeledSeqs:
    type: File
    outputBinding:
      glob: relabeledSeqs.fasta
  mergedMetadata:
    type: File
    outputBinding:
      glob: mergedMetadata.ttl
