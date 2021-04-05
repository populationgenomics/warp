version 1.0

import "WGSFromFastq.wdl" as WGSFromFastq
import "../../../../../../tasks/broad/BamProcessing.wdl" as Processing
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

workflow WGSFromBam {

  String pipeline_version = "2.1.0"

  input {
    File? input_cram
    File? input_bam
    File? output_map

    String sample_name
    String base_file_name
    String final_gvcf_base_name

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File wgs_coverage_interval_list

    Boolean to_cram = false
    Boolean check_contamination = true
    Boolean check_fingerprints = false
    Boolean provide_bam_output = false
    Boolean validate_gvcf = true
  }

  Int compression_level = 2
  Float input_size = size(input_bam, "GiB")
  # SortSam currently takes too long for preemptibles if the input data is too large
  Float gb_size_cutoff_for_preemptibles = 110.0
  Boolean input_too_large_for_preemptibles_1 = input_size > gb_size_cutoff_for_preemptibles

  call Processing.SortSamByName {
    input:
      input_bam = select_first([input_bam, input_cram]),
      output_bam_basename = base_file_name + ".sorted",
      compression_level = compression_level,
      preemptible_tries = if input_too_large_for_preemptibles_1 then 0 else papi_settings.agg_preemptible_tries,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
  }

  call BamToFastq {
    input:
      input_bam = SortSamByName.output_bam,
      output_fq_basename = base_file_name,
      disk_size = ceil(size(SortSamByName.output_bam, "GiB") * 6) + 20
  }

  SampleAndFastqs sample_and_fastqs = object {
     sample_name: sample_name,
     base_file_name: base_file_name,
     fastqs: [[BamToFastq.output_fq1, BamToFastq.output_fq2]],
     final_gvcf_base_name: final_gvcf_base_name,
   }

  call WGSFromFastq.WGSFromFastq {
    input:
      sample_and_fastqs = sample_and_fastqs,
      references = references,
      scatter_settings = scatter_settings,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      papi_settings = papi_settings,
      wgs_coverage_interval_list = wgs_coverage_interval_list,

      to_cram = to_cram,
      check_contamination = check_contamination,
      check_fingerprints = check_fingerprints,
      validate_gvcf = validate_gvcf,
      provide_bam_output = provide_bam_output
  }

  output {
    File? cross_check_fingerprints_metrics = WGSFromFastq.cross_check_fingerprints_metrics

    File? selfSM = WGSFromFastq.selfSM
    Float? contamination = WGSFromFastq.contamination

    File calculate_read_group_checksum_md5 = WGSFromFastq.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = WGSFromFastq.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = WGSFromFastq.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = WGSFromFastq.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = WGSFromFastq.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = WGSFromFastq.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = WGSFromFastq.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = WGSFromFastq.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = WGSFromFastq.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = WGSFromFastq.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = WGSFromFastq.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = WGSFromFastq.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = WGSFromFastq.agg_quality_distribution_metrics

    File? fingerprint_summary_metrics = WGSFromFastq.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = WGSFromFastq.fingerprint_detail_metrics

    File wgs_metrics = WGSFromFastq.wgs_metrics
    File raw_wgs_metrics = WGSFromFastq.raw_wgs_metrics

    File? duplicate_metrics = WGSFromFastq.duplicate_metrics
    File? output_bqsr_reports = WGSFromFastq.output_bqsr_reports

    File? gvcf_summary_metrics = WGSFromFastq.gvcf_summary_metrics
    File? gvcf_detail_metrics = WGSFromFastq.gvcf_detail_metrics

    File? output_bam = WGSFromFastq.output_bam
    File? output_bam_index = WGSFromFastq.output_bam_index

    File? output_cram = WGSFromFastq.output_cram
    File? output_cram_index = WGSFromFastq.output_cram_index
    File? output_cram_md5 = WGSFromFastq.output_cram_md5

    File? validate_cram_file_report = WGSFromFastq.validate_cram_file_report

    File output_vcf = WGSFromFastq.output_vcf
    File output_vcf_index = WGSFromFastq.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}


task BamToFastq {
  input {
    File input_bam
    String output_fq_basename
    Int disk_size
    File? ref_fasta
    File? ref_fasta_index
  }
  
  String ref_arg = if defined(ref_fasta) then "--reference ~{ref_fasta}" else ""
  command <<<
    samtools fastq ~{input_bam} \
    -1 ~{output_fq_basename}.1.fq \
    -2 ~{output_fq_basename}.2.fq \
    ~{ref_arg} \
    -1 ~{output_fq_basename}_R1.fastq.gz \
    -2 ~{output_fq_basename}_R2.fastq.gz \
    -0 /dev/null -s /dev/null
  >>>

  runtime {
    docker: "biocontainers/samtools:1.3.1"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4 GiB"
    preemptible: 3
  }

  output {
    File output_fq1 = "~{output_fq_basename}_R1.fastq.gz"
    File output_fq2 = "~{output_fq_basename}_R2.fastq.gz"
  }
}
