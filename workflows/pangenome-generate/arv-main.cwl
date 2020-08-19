cwlVersion: v1.1
class: Workflow
requirements:
  SubworkflowFeatureRequirement: {}
inputs:
  src_project: string
  metadataSchema: File
  exclude: File?
outputs:
  odgiGraph:
    type: File
    outputSource: pangenome-generate/odgiGraph
  odgiPNG:
    type: File
    outputSource: pangenome-generate/odgiPNG
  spoaGFA:
    type: File
    outputSource: pangenome-generate/spoaGFA
  odgiRDF:
    type: File
    outputSource: pangenome-generate/odgiRDF
  readsMergeDedup:
    type: File
    outputSource: pangenome-generate/readsMergeDedup
  mergedMetadata:
    type: File
    outputSource: pangenome-generate/mergedMetadata
  indexed_paths:
    type: File
    outputSource: pangenome-generate/indexed_paths
  colinear_components:
    type: Directory
    outputSource: pangenome-generate/colinear_components
steps:
  collect-seqs:
    run: collect-seqs.cwl
    in:
      src_project: src_project
      schema: metadataSchema
      exclude: exclude
    out: [relabeledSeqs, mergedMetadata]
  pangenome-generate:
    run: pangenome-generate_spoa.cwl
    in:
      seqs: collect-seqs/relabeledSeqs
      metadata: collect-seqs/mergedMetadata
      exclude: exclude
    out: [odgiGraph, odgiPNG, spoaGFA, odgiRDF, readsMergeDedup, mergedMetadata, indexed_paths, colinear_components]
