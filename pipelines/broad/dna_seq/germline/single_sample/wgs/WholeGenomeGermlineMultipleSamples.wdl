version 1.0

import "WholeGenomeGermlineSingleSample.wdl" as WholeGenomeGermlineSingleSample
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"


workflow WholeGenomeGermlineMultipleSamples {

  String pipeline_version = "1.0.0"

  input {
    File sample_map

    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings
    
    Boolean realign = true
    Boolean to_cram = true
    Boolean check_contamination = true
    Boolean check_fingerprints = false

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File wgs_coverage_interval_list

    Boolean output_alignment_file = true
    Boolean use_gatk3_haplotype_caller = false
    Boolean validate_gvcf = true

    String? subset_region
  }

  Array[Array[String]] sample_map_lines = read_tsv(sample_map)
  Int num_samples = length(sample_map_lines)  

  scatter (idx in range(num_samples)) {
    String sample_name = sample_map_lines[idx][0]
    File bam_or_cram_or_fastq1 = sample_map_lines[idx][1]
    File bai_or_crai_or_fastq2 = sample_map_lines[idx][2]
      
    Input inp = object {
       sample_name: sample_name,
       base_file_name: sample_name,
       bam_or_cram_or_fastq1: bam_or_cram_or_fastq1,
       bai_or_crai_or_fastq2: bai_or_crai_or_fastq2,
       final_gvcf_base_name: sample_name
    }  
  
    call WholeGenomeGermlineSingleSample.WholeGenomeGermlineSingleSample {
      input:
        inp = inp,
        references = references,
        scatter_settings = scatter_settings,
        papi_settings = papi_settings,

        realign = realign,
        to_cram = to_cram,
        check_contamination = check_contamination,
        check_fingerprints = check_fingerprints,

        fingerprint_genotypes_file = fingerprint_genotypes_file,
        fingerprint_genotypes_index = fingerprint_genotypes_index,

        wgs_coverage_interval_list = wgs_coverage_interval_list,

        output_alignment_file = output_alignment_file,
        use_gatk3_haplotype_caller = use_gatk3_haplotype_caller,
        validate_gvcf = validate_gvcf,
  
        subset_region = subset_region
    }
  }

  output {
    Array[File?] cross_check_fingerprints_metrics = select_all(WholeGenomeGermlineSingleSample.cross_check_fingerprints_metrics)

    Array[File?] selfSM = select_all(WholeGenomeGermlineSingleSample.selfSM)
    Array[Float?] contamination = select_all(WholeGenomeGermlineSingleSample.contamination)

    Array[File] agg_alignment_summary_metrics = select_all(WholeGenomeGermlineSingleSample.agg_alignment_summary_metrics)
    Array[File] agg_bait_bias_detail_metrics = select_all(WholeGenomeGermlineSingleSample.agg_bait_bias_detail_metrics)
    Array[File] agg_bait_bias_summary_metrics = select_all(WholeGenomeGermlineSingleSample.agg_bait_bias_summary_metrics)
    Array[File] agg_gc_bias_detail_metrics = select_all(WholeGenomeGermlineSingleSample.agg_gc_bias_detail_metrics)
    Array[File] agg_gc_bias_pdf = select_all(WholeGenomeGermlineSingleSample.agg_gc_bias_pdf)
    Array[File] agg_gc_bias_summary_metrics = select_all(WholeGenomeGermlineSingleSample.agg_gc_bias_summary_metrics)
    Array[File] agg_insert_size_histogram_pdf = select_all(WholeGenomeGermlineSingleSample.agg_insert_size_histogram_pdf)
    Array[File] agg_insert_size_metrics = select_all(WholeGenomeGermlineSingleSample.agg_insert_size_metrics)
    Array[File] agg_pre_adapter_detail_metrics = select_all(WholeGenomeGermlineSingleSample.agg_pre_adapter_detail_metrics)
    Array[File] agg_pre_adapter_summary_metrics = select_all(WholeGenomeGermlineSingleSample.agg_pre_adapter_summary_metrics)
    Array[File] agg_quality_distribution_pdf = select_all(WholeGenomeGermlineSingleSample.agg_quality_distribution_pdf)
    Array[File] agg_quality_distribution_metrics = select_all(WholeGenomeGermlineSingleSample.agg_quality_distribution_metrics)
    Array[File] agg_error_summary_metrics = select_all(WholeGenomeGermlineSingleSample.agg_error_summary_metrics)

    Array[File?] fingerprint_summary_metrics = select_all(WholeGenomeGermlineSingleSample.fingerprint_summary_metrics)
    Array[File?] fingerprint_detail_metrics = select_all(WholeGenomeGermlineSingleSample.fingerprint_detail_metrics)

    Array[File] wgs_metrics = select_all(WholeGenomeGermlineSingleSample.wgs_metrics)
    Array[File] raw_wgs_metrics = select_all(WholeGenomeGermlineSingleSample.raw_wgs_metrics)

    Array[File?] duplicate_metrics = select_all(WholeGenomeGermlineSingleSample.duplicate_metrics)

    Array[File?] gvcf_summary_metrics = select_all(WholeGenomeGermlineSingleSample.gvcf_summary_metrics)
    Array[File?] gvcf_detail_metrics = select_all(WholeGenomeGermlineSingleSample.gvcf_detail_metrics)

    Array[File?] alignment_file = select_all(WholeGenomeGermlineSingleSample.alignment_file)
    Array[File?] alignment_indx = select_all(WholeGenomeGermlineSingleSample.alignment_indx)

    Array[File] output_vcf = select_all(WholeGenomeGermlineSingleSample.output_vcf)
    Array[File] output_vcf_index = select_all(WholeGenomeGermlineSingleSample.output_vcf_index)
  }
  meta {
    allowNestedInputs: true
  }
}
