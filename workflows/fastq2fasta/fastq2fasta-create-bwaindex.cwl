cwlVersion: v1.1
class: Workflow
requirements:
  SubworkflowFeatureRequirement: {}
hints:
  ResourceRequirement:
    ramMin: 3000

inputs:
  ref_fasta:
    type: File
  fastq_forward:
    type: File
  fastq_reverse:
    type: File?
  threads:
    type: int
    default: 4

outputs:
  out_fasta:
    type: File
    outputSource: fastq2fasta/out_fasta

steps:
  bwa-index:
    in: {input_fasta: ref_fasta}
    out: [indexed_fasta]
    run: bwa-index.cwl
  samtools-faidx:
    in: {input_fasta: bwa-index/indexed_fasta}
    out: [indexed_fasta]
    run: samtools-faidx.cwl
  fastq2fasta:
    in:
      fastq_forward: fastq_forward
      fastq_reverse: fastq_reverse
      ref_fasta: samtools-faidx/indexed_fasta
      threads: threads
    out: [out_fasta]
    run: fastq2fasta.cwl
