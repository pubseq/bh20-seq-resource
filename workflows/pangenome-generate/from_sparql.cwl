cwlVersion: v1.1
class: CommandLineTool
requirements:
  DockerRequirement:
    dockerFile: |
      FROM debian:10
      RUN apt-get update && apt-get -yq --no-install-recommends install samtools python3-rdflib
    dockerImageId: rdflib-and-samtools
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
