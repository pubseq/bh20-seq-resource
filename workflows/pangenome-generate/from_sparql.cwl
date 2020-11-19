cwlVersion: v1.1
class: CommandLineTool
$namespaces:
  arv: "http://arvados.org/cwl#"
requirements:
  DockerRequirement:
    dockerFile: |
      FROM debian:10
      RUN apt-get update && apt-get -yq --no-install-recommends install samtools python3-rdflib
    dockerImageId: rdflib-and-samtools
  ResourceRequirement:
    ramMin: 768
  arv:RuntimeConstraints:
    keep_cache: 2048
    outputDirType: keep_output_dir
inputs:
  script:
    type: File
    default:
      class: File
      location: from_sparql.py
  metadata: File
  fasta:
    type: File
    secondaryFiles: [.fai]
  query: string
stdout: selected.fasta
outputs:
  selected: stdout
arguments: [python3, $(inputs.script), $(inputs.metadata), $(inputs.fasta), $(inputs.query)]
