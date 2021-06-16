version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human whole-genome data.
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in a BAM or CRAM format
## - GVCF output names must end in ".g.vcf.gz"
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "../../../../../../tasks/broad/UnmappedBamToAlignedBam.wdl" as ToBam
import "../../../../../../tasks/broad/AggregatedBamQC.wdl" as AggregatedQC
import "../../../../../../tasks/broad/Qc.wdl" as QC
import "../../../../../../tasks/broad/BamToCram.wdl" as ToCram
import "../../../../../../pipelines/broad/dna_seq/germline/variant_calling/VariantCalling.wdl" as ToGvcf
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow WholeGenomeGermlineSingleSample {

  String pipeline_version = "2.3.3"

  input {
    SampleAndBam sample_and_bam
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

  # Not overridable:
  Int read_length = 250
  Float lod_threshold = -20.0
  String cross_check_fingerprints_by = "READGROUP"
  String final_gvcf_base_name = select_first([sample_and_bam.final_gvcf_base_name, sample_and_bam.base_file_name])
  
  if (realign) {
    call ToBam.UnmappedBamToAlignedBam {
      input:
        sample_and_bam = sample_and_bam,
        references     = references,
        papi_settings  = papi_settings,
  
        check_contamination = check_contamination,
        check_fingerprints  = check_fingerprints,
  
        contamination_sites_ud  = references.contamination_sites_ud,
        contamination_sites_bed = references.contamination_sites_bed,
        contamination_sites_mu  = references.contamination_sites_mu,
  
        cross_check_fingerprints_by = cross_check_fingerprints_by,
        haplotype_database_file     = references.haplotype_database_file,
        lod_threshold               = lod_threshold,
      
        to_cram = to_cram,
        subset_region = subset_region,
    }
  }

  File mapped_file = select_first([UnmappedBamToAlignedBam.output_bam, sample_and_bam.bam_or_cram])
  File mapped_index = select_first([UnmappedBamToAlignedBam.output_index, sample_and_bam.bai_or_crai])

  call AggregatedQC.AggregatedBamQC {
    input:
      base_recalibrated_bam = mapped_file,
      base_recalibrated_bam_index = mapped_index,
      base_name = sample_and_bam.base_file_name,
      sample_name = sample_and_bam.sample_name,
      haplotype_database_file = references.haplotype_database_file,
      references = references,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      papi_settings = papi_settings
  }

  # QC the sample WGS metrics (stringent thresholds)
  call QC.CollectWgsMetrics as CollectWgsMetrics {
    input:
      input_bam = mapped_file,
      input_bam_index = mapped_index,
      metrics_filename = sample_and_bam.base_file_name + ".wgs_metrics",
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      read_length = read_length,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  # QC the sample raw WGS metrics (common thresholds)
  call QC.CollectRawWgsMetrics as CollectRawWgsMetrics {
    input:
      input_bam = mapped_file,
      input_bam_index = mapped_index,
      metrics_filename = sample_and_bam.base_file_name + ".raw_wgs_metrics",
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      read_length = read_length,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  call ToGvcf.VariantCalling as BamToGvcf {
    input:
      calling_interval_list = references.calling_interval_list,
      evaluation_interval_list = references.evaluation_interval_list,
      haplotype_scatter_count = scatter_settings.haplotype_scatter_count,
      break_bands_at_multiples_of = scatter_settings.break_bands_at_multiples_of,
      contamination = UnmappedBamToAlignedBam.contamination,
      input_bam = mapped_file,
      input_bam_index = mapped_index,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      ref_dict = references.reference_fasta.ref_dict,
      dbsnp_vcf = references.dbsnp_vcf,
      dbsnp_vcf_index = references.dbsnp_vcf_index,
      base_file_name = sample_and_bam.base_file_name,
      final_vcf_base_name = final_gvcf_base_name,
      agg_preemptible_tries = papi_settings.agg_preemptible_tries,
      use_gatk3_haplotype_caller = use_gatk3_haplotype_caller,
      validate_gvcf = validate_gvcf
  }

  if (output_alignment_file) {
    File alignment_file_to_output = mapped_file
    File alignment_index_to_output = mapped_index
  }

  # Outputs that will be retained when execution is complete
  output {
    File? cross_check_fingerprints_metrics = UnmappedBamToAlignedBam.cross_check_fingerprints_metrics

    File? selfSM = UnmappedBamToAlignedBam.selfSM
    Float? contamination = UnmappedBamToAlignedBam.contamination

    File calculate_read_group_checksum_md5 = AggregatedBamQC.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = AggregatedBamQC.agg_alignment_summary_metrics
    File agg_bait_bias_detail_metrics = AggregatedBamQC.agg_bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = AggregatedBamQC.agg_bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = AggregatedBamQC.agg_gc_bias_detail_metrics
    File agg_gc_bias_pdf = AggregatedBamQC.agg_gc_bias_pdf
    File agg_gc_bias_summary_metrics = AggregatedBamQC.agg_gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = AggregatedBamQC.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = AggregatedBamQC.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = AggregatedBamQC.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = AggregatedBamQC.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = AggregatedBamQC.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = AggregatedBamQC.agg_quality_distribution_metrics
    File agg_error_summary_metrics = AggregatedBamQC.agg_error_summary_metrics

    File? fingerprint_summary_metrics = AggregatedBamQC.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = AggregatedBamQC.fingerprint_detail_metrics

    File wgs_metrics = CollectWgsMetrics.metrics
    File raw_wgs_metrics = CollectRawWgsMetrics.metrics

    File? duplicate_metrics = UnmappedBamToAlignedBam.duplicate_metrics

    File? gvcf_summary_metrics = BamToGvcf.vcf_summary_metrics
    File? gvcf_detail_metrics = BamToGvcf.vcf_detail_metrics

    File? alignment_file = alignment_file_to_output
    File? alignment_index = alignment_index_to_output

    File output_vcf = BamToGvcf.output_vcf
    File output_vcf_index = BamToGvcf.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
