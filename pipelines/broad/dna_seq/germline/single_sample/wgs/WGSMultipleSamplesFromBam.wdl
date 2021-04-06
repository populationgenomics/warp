version 1.0

import "WGSFromBam.wdl" as WGSFromBam


workflow WGSMultipleSamplesFromBam {

  String pipeline_version = "1.0.0"

  input {
    File sample_map

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings
    Boolean to_cram = false
    Boolean check_contamination = true
    Boolean check_fingerprints = false
    Boolean do_bqsr = false

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File wgs_coverage_interval_list

    Boolean provide_bam_output = false
    Boolean use_gatk3_haplotype_caller = true
    Boolean validate_gvcf = true
  }

  Array[Array[String]] sample_map_lines = read_tsv(sample_map)
  Int num_samples = length(sample_map_lines)  

  scatter (idx in range(num_samples)) {
    String sample_name = sample_map_lines[idx][0]
    File input_bam = sample_map_lines[idx][1]
    
    call WGSFromBam.WGSFromBam as WGSFromBam {
      input:
        input_bam = input_bam,
        sample_name = sample_name,
        base_file_name = sample_name,
        final_gvcf_base_name = sample_name,
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
  }

  output {
    Array[File?] cross_check_fingerprints_metrics = select_all(WGSFromBam.cross_check_fingerprints_metrics)

    Array[File?] selfSM = select_all(WGSFromBam.selfSM)
    Array[Float?] contamination = select_all(WGSFromBam.contamination)

    Array[File] calculate_read_group_checksum_md5 = select_all(WGSFromBam.calculate_read_group_checksum_md5)

    Array[File] agg_alignment_summary_metrics = select_all(WGSFromBam.agg_alignment_summary_metrics)
    Array[File] agg_bait_bias_detail_metrics = select_all(WGSFromBam.agg_bait_bias_detail_metrics)
    Array[File] agg_bait_bias_summary_metrics = select_all(WGSFromBam.agg_bait_bias_summary_metrics)
    Array[File] agg_gc_bias_detail_metrics = select_all(WGSFromBam.agg_gc_bias_detail_metrics)
    Array[File] agg_gc_bias_pdf = select_all(WGSFromBam.agg_gc_bias_pdf)
    Array[File] agg_gc_bias_summary_metrics = select_all(WGSFromBam.agg_gc_bias_summary_metrics)
    Array[File] agg_insert_size_histogram_pdf = select_all(WGSFromBam.agg_insert_size_histogram_pdf)
    Array[File] agg_insert_size_metrics = select_all(WGSFromBam.agg_insert_size_metrics)
    Array[File] agg_pre_adapter_detail_metrics = select_all(WGSFromBam.agg_pre_adapter_detail_metrics)
    Array[File] agg_pre_adapter_summary_metrics = select_all(WGSFromBam.agg_pre_adapter_summary_metrics)
    Array[File] agg_quality_distribution_pdf = select_all(WGSFromBam.agg_quality_distribution_pdf)
    Array[File] agg_quality_distribution_metrics = select_all(WGSFromBam.agg_quality_distribution_metrics)

    Array[File?] fingerprint_summary_metrics = select_all(WGSFromBam.fingerprint_summary_metrics)
    Array[File?] fingerprint_detail_metrics = select_all(WGSFromBam.fingerprint_detail_metrics)

    Array[File] wgs_metrics = select_all(WGSFromBam.wgs_metrics)
    Array[File] raw_wgs_metrics = select_all(WGSFromBam.raw_wgs_metrics)

    Array[File?] duplicate_metrics = select_all(WGSFromBam.duplicate_metrics)
    Array[File?] output_bqsr_reports = select_all(WGSFromBam.output_bqsr_reports)

    Array[File?] gvcf_summary_metrics = select_all(WGSFromBam.gvcf_summary_metrics)
    Array[File?] gvcf_detail_metrics = select_all(WGSFromBam.gvcf_detail_metrics)

    Array[File?] output_bam = select_all(WGSFromBam.output_bam)
    Array[File?] output_bam_index = select_all(WGSFromBam.output_bam_index)

    Array[File?] output_cram = select_all(WGSFromBam.output_cram)
    Array[File?] output_cram_index = select_all(WGSFromBam.output_cram_index)
    Array[File?] output_cram_md5 = select_all(WGSFromBam.output_cram_md5)

    Array[File?] validate_cram_file_report = select_all(WGSFromBam.validate_cram_file_report)

    Array[File] output_vcf = select_all(WGSFromBam.output_vcf)
    Array[File] output_vcf_index = select_all(WGSFromBam.output_vcf_index)
  }
  meta {
    allowNestedInputs: true
  }
}
