#!/usr/bin/env cwl-runner
cwlVersion: v1.1
class: Workflow
requirements:
  ScatterFeatureRequirement: {}
  StepInputExpressionRequirement: {}
inputs:
  seqs: File
  metadata: File
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
#  odgiPNG:
#    type: File
#    outputSource: vizGraph/graph_image
  spoaGFA:
    type: File
    outputSource: induceGraph/spoaGFA
  odgiRDF:
    type: File
    outputSource: odgi2rdf/rdf
  readsMergeDedupSortedByQualAndLen:
    type: File
    outputSource: dedup_and_sort_by_quality_and_len/reads_dedupped_sorted_by_quality_and_len
  mergedMetadata:
    type: File
    outputSource: dups2metadata/merged
#  indexed_paths:
#    type: File
#    outputSource: index_paths/indexed_paths
#  colinear_components:
#    type: Directory
#    outputSource: segment_components/colinear_components
steps:
  dedup_and_sort_by_quality_and_len:
    in: {reads: seqs}
    out: [sortedReadsFA, dups]
    run: sort_fasta_by_quality_and_len.cwl
  induceGraph:
    in:
      readsFA: dedup_and_sort_by_quality_and_len/reads_dedupped_sorted_by_quality_and_len
    out: [spoaGFA]
    run: spoa.cwl
  buildGraph:
    in: {inputGFA: induceGraph/spoaGFA}
    out: [odgiGraph]
    run: odgi-build-from-spoa-gfa.cwl
  # vizGraph:
  #   in:
  #     sparse_graph_index: buildGraph/odgiGraph
  #     width:
  #       default: 50000
  #     height:
  #       default: 500
  #     path_per_row:
  #       default: true
  #     path_height:
  #       default: 4
  #   out: [graph_image]
  #   requirements:
  #     ResourceRequirement:
  #       ramMin: $(15 * 1024)
  #       outdirMin: 10
  #   run: ../tools/odgi/odgi_viz.cwl
  odgi2rdf:
    in: {odgi: buildGraph/odgiGraph}
    out: [rdf]
    run: odgi_to_rdf.cwl
  dups2metadata:
    in:
      metadata: metadata
      dups: dedup_and_sort_by_quality_and_len/dups
    out: [merged]
    run: dups2metadata.cwl
  # bin_paths:
  #   requirements:
  #     ResourceRequirement:
  #       ramMin: 3000
  #       outdirMin: 10
  #   run: ../tools/odgi/odgi_bin.cwl
  #   in:
  #     sparse_graph_index: buildGraph/odgiGraph
  #     bin_width: bin_widths
  #   scatter: bin_width
  #   out: [ bins, pangenome_sequence ]
  # index_paths:
  #   label: Create path index
  #   requirements:
  #     ResourceRequirement:
  #       ramMin: 3000
  #       outdirMin: 10
  #   run: ../tools/odgi/odgi_pathindex.cwl
  #   in:
  #     sparse_graph_index: buildGraph/odgiGraph
  #   out: [ indexed_paths ]
  # segment_components:
  #   label: Run component segmentation
  #   run: ../tools/graph-genome-segmentation/component_segmentation.cwl
  #   in:
  #     bins: bin_paths/bins
  #     cells_per_file: cells_per_file
  #     pangenome_sequence:
  #       source: bin_paths/pangenome_sequence
  #       valueFrom: $(self[0])
  #       # the bin_paths step is scattered over the bin_width array, but always using the same sparse_graph_index
  #       # the pangenome_sequence that is extracted is exactly the same for the same sparse_graph_index
  #       # regardless of bin_width, so we take the first pangenome_sequence as input for this step
  #   out: [ colinear_components ]
