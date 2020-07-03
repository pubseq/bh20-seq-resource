cwlVersion: v1.1
class: Workflow
inputs:
  existing_metadata_from_nuccore:
    type: Directory?
outputs: []
requirements:
  ScatterFeatureRequirement: {}
steps:
  fetch_from_genbank:
    in:
      existing_metadata_from_nuccore: existing_metadata_from_nuccore
    out: [fasta_and_yaml, metadata_from_nuccore, accessions]
    run: fetch_from_genbank.cwl
  split_into_arrays:
    in:
      dir: fetch_from_genbank/fasta_and_yaml
    out: [fasta, metadata]
    run: split_into_arrays.cwl
  upload:
    in:
      fasta: split_into_arrays/fasta
      metadata: split_into_arrays/metadata
    out: []
    scatter: [fasta, metadata]
    scatterMethod: dotproduct
    run: upload.cwl
