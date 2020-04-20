cwlVersion: v1.1
class: Workflow
requirements:
  SubworkflowFeatureRequirement: {}
hints:
  ResourceRequirement:
    ramMin: 3000

inputs:
  fastq_forward: File
  fastq_reverse: File?
  ref_fasta:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
      - .fai
  threads:
    type: int
    default: 4
  metadata: File?

outputs:
  out_fasta:
    type: File
    outputSource: bam2fasta/out_fasta
  out_metadata:
    type: File?
    outputSource: metadata

steps:
  bwa-mem:
    in:
      threads: threads
      fastq_forward: fastq_forward
      fastq_reverse: fastq_reverse
      index_base: ref_fasta
    out: [output]
    run: bwa-mem.cwl
  samtools-view:
    in:
      threads: threads
      input_file: bwa-mem/output
    out: [bam]
    run: samtools-view.cwl
  samtools-sort:
    in:
      input_bamfile: samtools-view/bam
      threads: threads
    out: [sorted_bam]
    run: samtools-sort.cwl
  bam2fasta:
    in:
      bam: samtools-sort/sorted_bam
      fasta: ref_fasta
      threads: threads
    out: [out_fasta]
    run: bam2fasta.cwl
