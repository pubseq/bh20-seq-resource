#!/usr/bin/env cwl-runner

cwlVersion: v1.1
class: Workflow

#############################################

inputs:

  # align
  threads: int
  pggb_wfmash: boolean
  pggb_fasta: File
  pggb_mash_k_mer: int
  pggb_map_percent_identity: int
  pggb_num_secondary_mappings: int
  pggb_segment_length: int
  pggb_output_dir: string

  # extract coverage vector
  odgi_paths: string
  odgi_graph: File
  haplotypes: boolean
  threads: int

  # remove consensus paths
  consensus_regex: string
  coverage_tsv: File

  # Get metadata
  main_py_script: File
  metadata: string
  coverage_matrix: File
  coverage_matrix_with_metadata: string

  # Generate newick tree
  main_py_script: File
  newick: string
  newick_dimensions: int
  newick_coverage_matrix: File  
  newick_metadata: string
  newick_tree: string

  # Genenrate augur JSON file
  nextstrain_bash_script: File
  newick_tree_2: File
  metadata_newick: File
  dataDir: Directory


#############################################

outputs:
  augur_json:
    type: File
    outputSource: augur/newick_json

#############################################

steps:
  align:
    run: align.cwl
    in:
      threads: threads
      pggb_wfmash: pggb_wfmash
      pggb_fasta: pggb_fasta
      pggb_mash_k_mer: pggb_mash_k_mer
      pggb_map_percent_identity: pggb_map_percent_identity
      pggb_num_secondary_mappings: pggb_num_secondary_mappings
      pggb_segment_length: pggb_segment_length
      pggb_output_dir: pggb_output_dir
    out: [pggb_odgi_graph]

  odgi:
    run: coverage.cwl
    in:
      odgi_paths: odgi_paths
      odgi_graph: align/pggb_odgi_graph
      haplotypes: haplotypes
      threads: threads
    out: [coverage_matrix]

  awk:
    run: awk-coverage.cwl
    in:
      consensus_regex: consensus_regex
      coverage_tsv: odgi/coverage_matrix
    out: [awk_coverage_matrix]

  metadata:
    run: metadata.cwl
    in:
      main_py_script: main_py_script
      metadata: metadata
      coverage_matrix: awk/awk_coverage_matrix
      coverage_matrix_with_metadata: coverage_matrix_with_metadata
    out: [coverage_matrix_with_metadata_out]

  newick:
    run: newick.cwl
    in:
      main_py_script: main_py_script
      newick: newick
      newick_dimensions: newick_dimensions
      newick_coverage_matrix: metadata/coverage_matrix_with_metadata_out
      newick_metadata: newick_metadata
      newick_tree: newick_tree      
    out: [metadata_out, newick_tree_out]

  augur:
    run: augur.cwl
    in:
      nextstrain_bash_script: nextstrain_bash_script
      newick_tree_2: newick/newick_tree_out
      metadata_newick: newick/metadata_out
      dataDir: dataDir
      
    out: [newick_json]
