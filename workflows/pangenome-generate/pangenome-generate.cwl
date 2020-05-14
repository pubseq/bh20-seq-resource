#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
inputs:
  inputReads: File[]
  metadata: File[]
  metadataSchema: File
  subjects: string[]
  bin_widths:
    type: int[]
    default: [ 1, 4, 16, 64, 256, 1000, 4000, 16000]
    doc: width of each bin in basepairs along the graph vector
  cells_per_file:
    type: int
    default: 100
    doc: Cells per file on component_segmentation
outputs:
  odgiGraph:
    type: File
    outputSource: buildGraph/odgiGraph
  odgiPNG:
    type: File
    outputSource: vizGraph/graph_image
  seqwishGFA:
    type: File
    outputSource: induceGraph/seqwishGFA
  odgiRDF:
    type: File
    outputSource: odgi2rdf/rdf
  readsMergeDedup:
    type: File
    outputSource: dedup/reads_dedup
  mergedMetadata:
    type: File
    outputSource: mergeMetadata/merged
  indexed_paths:
    type: File
    outputSource: index_paths/indexed_paths
  colinear_components:
    type: Directory
    outputSource: segment_components/colinear_components
steps:
  relabel:
    in:
      readsFA: inputReads
      subjects: subjects
    out: [relabeledSeqs, originalLabels]
    run: relabel-seqs.cwl
  dedup:
    in: {reads: relabel/relabeledSeqs}
    out: [reads_dedup, dups]
    run: ../tools/seqkit/seqkit_rmdup.cwl
  overlapReads:
    in:
      target: dedup/reads_dedup
      query: dedup/reads_dedup
      outputCIGAR: {default: true}
      preset: {default: asm20}
      miniWinSize: {default: 1}
    out: [alignments]
    run: ../tools/minimap2/minimap2_paf.cwl
  induceGraph:
    in:
      readsFA: dedup/reads_dedup
      readsPAF: overlapReads/alignments
    out: [seqwishGFA]
    run: seqwish.cwl
  buildGraph:
    in: {inputGFA: induceGraph/seqwishGFA}
    out: [odgiGraph]
    run: odgi-build.cwl
  vizGraph:
    in:
      sparse_graph_index: buildGraph/odgiGraph
      width:
        default: 50000
      height:
        default: 500
      path_per_row:
        default: true
      path_height:
        default: 4
    out: [graph_image]
    run: ../tools/odgi/odgi_viz.cwl
  odgi2rdf:
    in: {odgi: buildGraph/odgiGraph}
    out: [rdf]
    run: odgi_to_rdf.cwl
  mergeMetadata:
    in:
      metadata: metadata
      metadataSchema: metadataSchema
      subjects: subjects
      dups: dedup/dups
      originalLabels: relabel/originalLabels
    out: [merged]
    run: merge-metadata.cwl
  bin_paths:
    run: ../tools/odgi/odgi_bin.cwl
    in:
      sparse_graph_index: buildGraph/odgiGraph
      bin_width: bin_widths
    scatter: bin_width
    out: [ bins, pangenome_sequence ]
  index_paths:
    label: Create path index
    run: ../tools/odgi/odgi_pathindex.cwl
    in:
      sparse_graph_index: buildGraph/odgiGraph
    out: [ indexed_paths ]
  segment_components:
    label: Run component segmentation
    run: ../tools/graph-genome-segmentation/component_segmentation.cwl
    in:
      bins: bin_paths/bins
      cells_per_file: cells_per_file
      pangenome_sequence:
        source: bin_paths/pangenome_sequence
        valueFrom: $(self[0])
        # the bin_paths step is scattered over the bin_width array, but always using the same sparse_graph_index
        # the pangenome_sequence that is extracted is exactly the same for the same sparse_graph_index
        # regardless of bin_width, so we take the first pangenome_sequence as input for this step
    out: [ colinear_components ]
