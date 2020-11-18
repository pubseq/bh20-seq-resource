cwlVersion: v1.1
class: Workflow
requirements:
  SubworkflowFeatureRequirement: {}
inputs:
  metadata: File
  fasta:
    type: File
    secondaryFiles: [.fai]
  query: string
outputs:
  odgiGraph:
    type: File
    outputSource: make-gfa/odgiGraph
  spoaGFA:
    type: File
    outputSource: make-gfa/spoaGFA
  readsMergeDedupSortedByQualAndLen:
    type: File
    outputSource: make-gfa/readsMergeDedupSortedByQualAndLen
  mergedMetadata:
    type: File
    outputSource: make-gfa/mergedMetadata
steps:
  get-subset:
    run: from_sparql.cwl
    in: {metadata: metadata, query: query, fasta: fasta}
    out: [selected]
  make-gfa:
    run: pangenome-generate_spoa.cwl
    in: {metadata: metadata, seqs: get-subset/selected}
    out: [odgiGraph, spoaGFA, readsMergeDedupSortedByQualAndLen, mergedMetadata]
