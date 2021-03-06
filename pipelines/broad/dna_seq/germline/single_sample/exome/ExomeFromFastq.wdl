version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data pre-processing and initial variant calling (GVCF
## generation) according to the GATK Best Practices (June 2016) for germline SNP and
## Indel discovery in human exome sequencing data.
##
## Requirements/expectations :
## - Human exome sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single sample (SM)
## - Input uBAM files must additionally comply with the following requirements:
## - - filenames all have the same suffix (we use ".unmapped.bam")
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
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

import "../../../../../../tasks/broad/FastqsToAlignedBam.wdl" as ToBam
import "../../../../../../tasks/broad/AggregatedBamQC.wdl" as AggregatedQC
import "../../../../../../tasks/broad/BamProcessing.wdl" as Processing
import "../../../../../../tasks/broad/Qc.wdl" as QC
import "../../../../../../tasks/broad/BamToCram.wdl" as ToCram
import "../../../../../../tasks/broad/VariantCalling.wdl" as ToGvcf
import "../../../../../../structs/dna_seq/DNASeqStructs.wdl"

# WORKFLOW DEFINITION
workflow ExomeFromFastq {


  String pipeline_version = "2.4.1"

  input {
    SampleAndFastqs sample_and_fastqs
    DNASeqSingleSampleReferences references
    VariantCallingScatterSettings scatter_settings
    PapiSettings papi_settings
    Boolean to_cram = false
    Boolean check_contamination = true
    Boolean check_fingerprints = false
    Boolean do_bqsr = false

    File? fingerprint_genotypes_file
    File? fingerprint_genotypes_index

    File target_interval_list

    Boolean provide_bam_output = false
    Boolean validate_gvcf = true
  }

  # Not overridable:
  Float lod_threshold = -10.0
  String cross_check_fingerprints_by = "READGROUP"
  String recalibrated_bam_basename = sample_and_fastqs.base_file_name + ".aligned.duplicates_marked.recalibrated"

  String final_gvcf_base_name = select_first([sample_and_fastqs.final_gvcf_base_name, sample_and_fastqs.base_file_name])

  call Processing.GenerateSubsettedContaminationResources {
    input:
        target_interval_list = target_interval_list,
        contamination_sites_bed = references.contamination_sites_bed,
        contamination_sites_mu = references.contamination_sites_mu,
        contamination_sites_ud = references.contamination_sites_ud,
        preemptible_tries = papi_settings.preemptible_tries
  }

  call ToBam.FastqsToAlignedBam {
    input:
      sample_and_fastqs = sample_and_fastqs,
      references     = references,
      papi_settings  = papi_settings,

      check_contamination = check_contamination,
      check_fingerprints  = check_fingerprints,
      do_bqsr             = do_bqsr,

      contamination_sites_ud = GenerateSubsettedContaminationResources.subsetted_contamination_ud,
      contamination_sites_bed = GenerateSubsettedContaminationResources.subsetted_contamination_bed,
      contamination_sites_mu = GenerateSubsettedContaminationResources.subsetted_contamination_mu,

      cross_check_fingerprints_by = cross_check_fingerprints_by,
      haplotype_database_file     = references.haplotype_database_file,
      lod_threshold               = lod_threshold,
      recalibrated_bam_basename   = recalibrated_bam_basename
  }

  File mapped_bam = FastqsToAlignedBam.output_bam
  File mapped_index = FastqsToAlignedBam.output_bam_index
  File? dup_metrics = FastqsToAlignedBam.duplicate_metrics

  call AggregatedQC.AggregatedBamQC {
    input:
      base_recalibrated_bam = mapped_bam,
      base_recalibrated_bam_index = mapped_index,
      base_name = sample_and_fastqs.base_file_name,
      sample_name = sample_and_fastqs.sample_name,
      recalibrated_bam_base_name = recalibrated_bam_basename,
      haplotype_database_file = references.haplotype_database_file,
      references = references,
      fingerprint_genotypes_file = fingerprint_genotypes_file,
      fingerprint_genotypes_index = fingerprint_genotypes_index,
      papi_settings = papi_settings
  }

  if (to_cram) {
    call ToCram.BamToCram as BamToCram {
      input:
        input_bam = mapped_bam,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        ref_dict = references.reference_fasta.ref_dict,
        duplication_metrics = dup_metrics,
        chimerism_metrics = AggregatedBamQC.agg_alignment_summary_metrics,
        base_file_name = sample_and_fastqs.base_file_name,
        agg_preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  call ToGvcf.VariantCalling as BamToGvcf {
    input:
      calling_interval_list = references.calling_interval_list,
      evaluation_interval_list = references.evaluation_interval_list,
      haplotype_scatter_count = scatter_settings.haplotype_scatter_count,
      break_bands_at_multiples_of = scatter_settings.break_bands_at_multiples_of,
      contamination = FastqsToAlignedBam.contamination,
      input_bam = mapped_bam,
      input_bai = mapped_index,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      ref_dict = references.reference_fasta.ref_dict,
      dbsnp_vcf = references.dbsnp_vcf,
      dbsnp_vcf_index = references.dbsnp_vcf_index,
      base_file_name = sample_and_fastqs.base_file_name,
      final_vcf_base_name = final_gvcf_base_name,
      agg_preemptible_tries = papi_settings.agg_preemptible_tries
  }

  call QC.CollectHsMetrics as CollectHsMetrics {
    input:
      input_bam = mapped_bam,
      input_bam_index = mapped_index,
      metrics_filename = sample_and_fastqs.base_file_name + ".hybrid_selection_metrics",
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      target_interval_list = target_interval_list,
      bait_interval_list = target_interval_list,
      preemptible_tries = papi_settings.agg_preemptible_tries
  }

  if (provide_bam_output) {
    File provided_output_bam = mapped_bam
    File provided_output_bam_index = mapped_index
  }

  # Outputs that will be retained when execution is complete
  output {
    File? cross_check_fingerprints_metrics = FastqsToAlignedBam.cross_check_fingerprints_metrics

    File? selfSM = FastqsToAlignedBam.selfSM
    Float? contamination = FastqsToAlignedBam.contamination

    File calculate_read_group_checksum_md5 = AggregatedBamQC.calculate_read_group_checksum_md5

    File agg_alignment_summary_metrics = AggregatedBamQC.agg_alignment_summary_metrics
    File agg_insert_size_histogram_pdf = AggregatedBamQC.agg_insert_size_histogram_pdf
    File agg_insert_size_metrics = AggregatedBamQC.agg_insert_size_metrics
    File agg_pre_adapter_detail_metrics = AggregatedBamQC.agg_pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = AggregatedBamQC.agg_pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = AggregatedBamQC.agg_quality_distribution_pdf
    File agg_quality_distribution_metrics = AggregatedBamQC.agg_quality_distribution_metrics
    File agg_error_summary_metrics = AggregatedBamQC.agg_error_summary_metrics

    File? fingerprint_summary_metrics = AggregatedBamQC.fingerprint_summary_metrics
    File? fingerprint_detail_metrics = AggregatedBamQC.fingerprint_detail_metrics

    File? duplicate_metrics = FastqsToAlignedBam.duplicate_metrics
    File? output_bqsr_reports = FastqsToAlignedBam.output_bqsr_reports

    File? gvcf_summary_metrics = BamToGvcf.vcf_summary_metrics
    File? gvcf_detail_metrics = BamToGvcf.vcf_detail_metrics

    File hybrid_selection_metrics = CollectHsMetrics.metrics

    File? output_bam = provided_output_bam
    File? output_bam_index = provided_output_bam_index

    File? output_cram = BamToCram.output_cram
    File? output_cram_index = BamToCram.output_cram_index
    File? output_cram_md5 = BamToCram.output_cram_md5

    File? validate_cram_file_report = BamToCram.validate_cram_file_report

    File output_vcf = BamToGvcf.output_vcf
    File output_vcf_index = BamToGvcf.output_vcf_index
  }
  meta {
    allowNestedInputs: true
  }
}
