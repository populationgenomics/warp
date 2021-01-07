version 1.0

import "ExomeFromFastq.wdl" as ExomeFromFastq
import "../../../../../../tasks/broad/BamProcessing.wdl" as Processing
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow ExomeFromBam {

  String pipeline_version = "2.1.0"

  input {
    File? input_cram
    File? input_bam
    File? output_map

    String sample_name
    String base_file_name
    String final_gvcf_base_name

    File? cram_ref_fasta
    File? cram_ref_fasta_index

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File target_interval_list

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
      preemptible_tries = if input_too_large_for_preemptibles_1 then 0 else papi_settings.agg_preemptible_tries
  }

  call BamToFastq {
    input:
      input_bam = SortSamByName.output_bam,
      output_fq_basename = base_file_name,
      disk_size = ceil(input_size * 6) + 20
  }

  SampleAndFastqs sample_and_fastqs = object {
     sample_name: sample_name,
     base_file_name: base_file_name,
     fastqs: [[BamToFastq.output_fq1, BamToFastq.output_fq2]],
     final_gvcf_base_name: final_gvcf_base_name,
  }

  call ExomeFromFastq.ExomeFromFastq {
    input:
      sample_and_fastqs = sample_and_fastqs,
      references = references,
      scatter_settings = scatter_settings,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      papi_settings = papi_settings,
      target_interval_list = target_interval_list,

      to_cram = to_cram,
      check_contamination = check_contamination,
      check_fingerprints = check_fingerprints,
      validate_gvcf = validate_gvcf,
      provide_bam_output = provide_bam_output
  }

  # Outputs that will be retained when execution is complete
  output {
#    Array[File] quality_yield_metrics = ExomeFromFastq.quality_yield_metrics
#
#    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = ExomeFromFastq.unsorted_read_group_base_distribution_by_cycle_pdf
#    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = ExomeFromFastq.unsorted_read_group_base_distribution_by_cycle_metrics
#    Array[File] unsorted_read_group_insert_size_histogram_pdf = ExomeFromFastq.unsorted_read_group_insert_size_histogram_pdf
#    Array[File] unsorted_read_group_insert_size_metrics = ExomeFromFastq.unsorted_read_group_insert_size_metrics
#    Array[File] unsorted_read_group_quality_by_cycle_pdf = ExomeFromFastq.unsorted_read_group_quality_by_cycle_pdf
#    Array[File] unsorted_read_group_quality_by_cycle_metrics = ExomeFromFastq.unsorted_read_group_quality_by_cycle_metrics
#    Array[File] unsorted_read_group_quality_distribution_pdf = ExomeFromFastq.unsorted_read_group_quality_distribution_pdf
#    Array[File] unsorted_read_group_quality_distribution_metrics = ExomeFromFastq.unsorted_read_group_quality_distribution_metrics
#
#    File read_group_alignment_summary_metrics = ExomeFromFastq.read_group_alignment_summary_metrics

    File? cross_check_fingerprints_metrics = ExomeFromFastq.cross_check_fingerprints_metrics

    File? selfSM = ExomeFromFastq.selfSM
    Float? contamination = ExomeFromFastq.contamination

    File calculate_read_group_checksum_md5 = ExomeFromFastq.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = ExomeFromFastq.agg_alignment_summary_metrics
#    File agg_bait_bias_detail_metrics = ExomeFromFastq.agg_bait_bias_detail_metrics
#    File agg_bait_bias_summary_metrics = ExomeFromFastq.agg_bait_bias_summary_metrics
    File agg_insert_size_histogram_pdf = ExomeFromFastq.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = ExomeFromFastq.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = ExomeFromFastq.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = ExomeFromFastq.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = ExomeFromFastq.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = ExomeFromFastq.agg_quality_distribution_metrics
    File agg_error_summary_metrics = ExomeFromFastq.agg_error_summary_metrics

    File? fingerprint_summary_metrics = ExomeFromFastq.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = ExomeFromFastq.fingerprint_detail_metrics

    File? duplicate_metrics = ExomeFromFastq.duplicate_metrics
    File? output_bqsr_reports = ExomeFromFastq.output_bqsr_reports

    File? gvcf_summary_metrics = ExomeFromFastq.gvcf_summary_metrics
    File? gvcf_detail_metrics = ExomeFromFastq.gvcf_detail_metrics

    File hybrid_selection_metrics = ExomeFromFastq.hybrid_selection_metrics

    File? output_bam = ExomeFromFastq.output_bam
    File? output_bam_index = ExomeFromFastq.output_bam_index

    File? output_cram = ExomeFromFastq.output_cram
    File? output_cram_index = ExomeFromFastq.output_cram_index
    File? output_cram_md5 = ExomeFromFastq.output_cram_md5

    File? validate_cram_file_report = ExomeFromFastq.validate_cram_file_report

    File output_vcf = ExomeFromFastq.output_vcf
    File output_vcf_index = ExomeFromFastq.output_vcf_index
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
  }

  command <<<
    samtools fastq ~{input_bam} \
    -1 ~{output_fq_basename}.1.fq \
    -2 ~{output_fq_basename}.2.fq \
    -0 /dev/null -s /dev/null
  >>>

  runtime {
    docker: "biocontainers/samtools:1.3.1"
    disks: "local-disk " + disk_size + " HDD"
    memory: "4 GiB"
    preemptible: 3
  }

  output {
    File output_fq1 = output_fq_basename + '.1.fq'
    File output_fq2 = output_fq_basename + '.2.fq'
  }
}
