version 1.0

## Copyright Broad Institute, 2018
##
## This WDL pipeline implements data processing according to the GATK Best Practices (June 2016)
## for human whole-genome and exome sequencing data.
##
## Runtime parameters are often optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

import "../../tasks/broad/Alignment.wdl" as Alignment
import "../../tasks/broad/SplitLargeReadGroup.wdl" as SplitRG
import "../../tasks/broad/Qc.wdl" as QC
import "../../tasks/broad/BamProcessing.wdl" as Processing
import "../../tasks/broad/Utilities.wdl" as Utils
import "../../structs/dna_seq/DNASeqStructs.wdl" as Structs

task BamToUnmappedBam {
  input {
    File input_bam
    String output_bam_filename
    Int disk_size

    
  }
  Int disk_size = ceil(size(input_bam, "GiB")) + 20
  disk_size = ceil(input_size * 3) + additional_disk

  command <<<

    java -Xms8G -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
      RevertSam \
      -I=~{input_bam} \
      -O=~{output_bam_filename} \
      --SANITIZE=true \
      --MAX_DISCARD_FRACTION=0.005 \
      --ATTRIBUTE_TO_CLEAR=XT \
      --ATTRIBUTE_TO_CLEAR=XN \
      --ATTRIBUTE_TO_CLEAR=AS \
      --ATTRIBUTE_TO_CLEAR=OC \
      --ATTRIBUTE_TO_CLEAR=OP \
      --ATTRIBUTE_TO_CLEAR=CO \
      --SORT_ORDER=queryname \
      --RESTORE_ORIGINAL_QUALITIES=true \
      --REMOVE_DUPLICATE_INFORMATION=true \
      --REMOVE_ALIGNMENT_INFORMATION=true \
      --VALIDATION_STRINGENCY=LENIENT

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.22.3"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.75 GiB"
    preemptible: 3
  }

  output {
    File output_bam = output_bam_filename
  }
}

# WORKFLOW DEFINITION
workflow UnmappedBamToAlignedBam {

  input {
    SampleAndBam sample_and_bam
    DNASeqSingleSampleReferences references
    PapiSettings papi_settings

    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    String cross_check_fingerprints_by
    File haplotype_database_file
    Float lod_threshold
    String recalibrated_bam_basename
    Boolean hard_clip_reads = false
    Boolean bin_base_qualities = true
    Boolean somatic = false
  }

  Float cutoff_for_large_rg_in_gb = 20.0

  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"

  Int compression_level = 2

  # Get the size of the standard reference files as well as the additional reference files needed for BWA

  String input_bam_basename = basename(sample_and_bam.mapped_bam, ".bam")

  call BamToUnmappedBam {
    input:
      input_bam = sample_and_bam.mapped_bam,
      output_bam_filename = input_bam_basename + ".unmapped.bam"
  }

  File unmapped_bam = BamToUnmappedBam.output_bam

  Float unmapped_bam_size = size(unmapped_bam, "GiB")

  # QC the unmapped BAM
  call QC.CollectQualityYieldMetrics as CollectQualityYieldMetrics {
    input:
      input_bam = unmapped_bam,
      metrics_filename = input_bam_basename + ".unmapped.quality_yield_metrics",
      preemptible_tries = papi_settings.preemptible_tries
  }

  if (unmapped_bam_size > cutoff_for_large_rg_in_gb) {
    # Split bam into multiple smaller bams,
    # map reads to reference and recombine into one bam
    call SplitRG.SplitLargeReadGroup as SplitRG {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = input_bam_basename + ".aligned.unsorted",
        reference_fasta = references.reference_fasta,
        compression_level = compression_level,
        preemptible_tries = papi_settings.preemptible_tries,
        hard_clip_reads = hard_clip_reads
    }
  }

  if (unmapped_bam_size <= cutoff_for_large_rg_in_gb) {
    # Map reads to reference
    call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = input_bam_basename + ".aligned.unsorted",
        reference_fasta = references.reference_fasta,
        compression_level = compression_level,
        preemptible_tries = papi_settings.preemptible_tries,
        hard_clip_reads = hard_clip_reads
    }
  }

  File output_aligned_bam = select_first([SamToFastqAndBwaMemAndMba.output_bam, SplitRG.aligned_bam])

  Float mapped_bam_size = size(output_aligned_bam, "GiB")

  # QC the aligned but unsorted readgroup BAM
  # no reference as the input here is unsorted, providing a reference would cause an error
  call QC.CollectUnsortedReadgroupBamQualityMetrics as CollectUnsortedReadgroupBamQualityMetrics {
    input:
      input_bam = output_aligned_bam,
      output_bam_prefix = input_bam_basename + ".readgroup",
      preemptible_tries = papi_settings.preemptible_tries
  }

  # Sum the read group bam sizes to approximate the aggregated bam size
  call Utils.SumFloats as SumFloats {
    input:
      sizes = [mapped_bam_size],
      preemptible_tries = papi_settings.preemptible_tries
  }

  # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
  Float gb_size_cutoff_for_preemptibles = 110.0
  Boolean data_too_large_for_preemptibles = SumFloats.total_size > gb_size_cutoff_for_preemptibles

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call Processing.MarkDuplicates as MarkDuplicates {
    input:
      input_bams = [output_aligned_bam],
      output_bam_basename = sample_and_bam.base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = sample_and_bam.base_file_name + ".duplicate_metrics",
      total_input_size = SumFloats.total_size,
      compression_level = compression_level,
      preemptible_tries = if data_too_large_for_preemptibles then 0 else papi_settings.agg_preemptible_tries
  }

  # Sort aggregated+deduped BAM file and fix tags
  call Processing.SortSam as SortSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = sample_and_bam.base_file_name + ".aligned.duplicate_marked.sorted",
      compression_level = compression_level,
      preemptible_tries = if data_too_large_for_preemptibles then 0 else papi_settings.agg_preemptible_tries
  }

  Float agg_bam_size = size(SortSampleBam.output_bam, "GiB")

  if (defined(haplotype_database_file)) {
    # Check identity of fingerprints across readgroups
    call QC.CrossCheckFingerprints as CrossCheckFingerprints {
      input:
        input_bams = [ SortSampleBam.output_bam ],
        input_bam_indexes = [SortSampleBam.output_bam_index],
        haplotype_database_file = haplotype_database_file,
        metrics_filename = sample_and_bam.base_file_name + ".crosscheck",
        total_input_size = agg_bam_size,
        lod_threshold = lod_threshold,
        cross_check_by = cross_check_fingerprints_by,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  # Create list of sequences for scatter-gather parallelization
  call Utils.CreateSequenceGroupingTSV as CreateSequenceGroupingTSV {
    input:
      ref_dict = references.reference_fasta.ref_dict,
      preemptible_tries = papi_settings.preemptible_tries
  }

  # Estimate level of cross-sample contamination
  call Processing.CheckContamination as CheckContamination {
    input:
      input_bam = SortSampleBam.output_bam,
      input_bam_index = SortSampleBam.output_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = references.reference_fasta.ref_fasta,
      ref_fasta_index = references.reference_fasta.ref_fasta_index,
      output_prefix = sample_and_bam.base_file_name,
      preemptible_tries = papi_settings.agg_preemptible_tries,
      contamination_underestimation_factor = 0.75
  }

  # Outputs that will be retained when execution is complete
  output {
    File quality_yield_metrics = CollectQualityYieldMetrics.quality_yield_metrics

    File unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf
    File unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics
    File unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf
    File unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics
    File unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf
    File unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics
    File unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf
    File unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics

    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination

    File duplicate_metrics = MarkDuplicates.duplicate_metrics

    File output_bam = SortSampleBam.output_bam
    File output_bam_index = SortSampleBam.output_bam_index
  }
  meta {
    allowNestedInputs: true
  }
}
