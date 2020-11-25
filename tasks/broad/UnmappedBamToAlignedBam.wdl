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
import "../../tasks/broad/Qc.wdl" as QC
import "../../tasks/broad/BamProcessing.wdl" as Processing
import "../../tasks/broad/Utilities.wdl" as Utils
import "../../structs/dna_seq/DNASeqStructs.wdl" as Structs

# WORKFLOW DEFINITION
workflow UnmappedBamToAlignedBam {

  input {
    SampleAndBam sample_and_bam
    DNASeqSingleSampleReferences references
    PapiSettings papi_settings

    Boolean hard_clip_reads = false
    Boolean bin_base_qualities = true
    Boolean somatic = false

    Boolean check_contamination = true
    Boolean check_fingerprints = false
    Boolean do_bqsr = false

    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    String cross_check_fingerprints_by
    File haplotype_database_file
    Float lod_threshold
    String recalibrated_bam_basename
  }

  Float cutoff_for_large_rg_in_gb = 20.0

  String sample_name = sample_and_bam.sample_name
  String rg_line = "@RG\\tID:" + sample_name + "\\tSM:" + sample_name
  String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y -R '" + rg_line + "' $bash_ref_fasta "

  Int compression_level = 2

  # Get the size of the standard reference files as well as the additional reference files needed for BWA

  String input_bam_basename = basename(sample_and_bam.mapped_bam, ".bam")

  Float input_size = size(sample_and_bam.mapped_bam, "GiB")
  # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
  Float gb_size_cutoff_for_preemptibles = 110.0
  Boolean input_too_large_for_preemptibles_1 = input_size > gb_size_cutoff_for_preemptibles

  call Processing.SortSamByName {
    input:
      input_bam = sample_and_bam.mapped_bam,
      output_bam_basename = input_bam_basename + ".sorted",
      compression_level = compression_level,
      preemptible_tries = if input_too_large_for_preemptibles_1 then 0 else papi_settings.agg_preemptible_tries
  }

  call BamToFastq {
    input:
      input_bam = SortSamByName.output_bam,
      output_fq_basename = input_bam_basename,
      disk_size = ceil(input_size * 6) + 20
  }

  File fq1 = BamToFastq.output_fq1
  File fq2 = BamToFastq.output_fq2
  call Alignment.SamToFastqAndBwaMemAndMba as SamToFastqAndBwaMemAndMba {
    input:
      input_fq1 = fq1,
      input_fq2 = fq2,
      bwa_commandline = bwa_commandline,
      output_bam_basename = input_bam_basename + ".aligned.unsorted",
      reference_fasta = references.reference_fasta,
      compression_level = compression_level,
      preemptible_tries = papi_settings.preemptible_tries,
      hard_clip_reads = hard_clip_reads
  }
#  }

  File output_aligned_bam = SamToFastqAndBwaMemAndMba.output_bam
  Float mapped_bam_size = size(output_aligned_bam, "GiB")

  # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
  Boolean input_too_large_for_preemptibles_2 = mapped_bam_size > gb_size_cutoff_for_preemptibles

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call Processing.MarkDuplicates as MarkDuplicates {
    input:
      input_bams = [output_aligned_bam],
      output_bam_basename = sample_and_bam.base_file_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = sample_and_bam.base_file_name + ".duplicate_metrics",
      total_input_size = mapped_bam_size,
      compression_level = compression_level,
      preemptible_tries = if input_too_large_for_preemptibles_2 then 0 else papi_settings.agg_preemptible_tries
  }

  # Sort aggregated+deduped BAM file and fix tags
  call Processing.SortSam as SortSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = sample_and_bam.base_file_name + ".aligned.duplicate_marked.sorted",
      compression_level = compression_level,
      preemptible_tries = if input_too_large_for_preemptibles_2 then 0 else papi_settings.agg_preemptible_tries
  }

  Float agg_bam_size = size(SortSampleBam.output_bam, "GiB")

  if (check_fingerprints && defined(haplotype_database_file)) {
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

  if (check_contamination) {
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
  }

  if (do_bqsr && defined(references.dbsnp_vcf)) {
    # We need disk to localize the sharded input and output due to the scatter for BQSR.
    # If we take the number we are scattering by and reduce by 3 we will have enough disk space
    # to account for the fact that the data is not split evenly.
    Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
    Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
    Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

    # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
    scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
      # Generate the recalibration model by interval
      call Processing.BaseRecalibrator as BaseRecalibrator {
        input:
          input_bam = SortSampleBam.output_bam,
          recalibration_report_filename = sample_and_bam.base_file_name + ".recal_data.csv",
          sequence_group_interval = subgroup,
          dbsnp_vcf = select_first([references.dbsnp_vcf]),
          dbsnp_vcf_index = select_first([references.dbsnp_vcf_index]),
          known_indels_sites_vcfs = references.known_indels_sites_vcfs,
          known_indels_sites_indices = references.known_indels_sites_indices,
          ref_dict = references.reference_fasta.ref_dict,
          ref_fasta = references.reference_fasta.ref_fasta,
          ref_fasta_index = references.reference_fasta.ref_fasta_index,
          bqsr_scatter = bqsr_divisor,
          preemptible_tries = papi_settings.agg_preemptible_tries
      }
    }

    # Merge the recalibration reports resulting from by-interval recalibration
    # The reports are always the same size
    call Processing.GatherBqsrReports as GatherBqsrReports {
      input:
        input_bqsr_reports = BaseRecalibrator.recalibration_report,
        output_report_filename = sample_and_bam.base_file_name + ".recal_data.csv",
        preemptible_tries = papi_settings.preemptible_tries
    }

    scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
      # Apply the recalibration model by interval
      call Processing.ApplyBQSR as ApplyBQSR {
        input:
          input_bam = SortSampleBam.output_bam,
          output_bam_basename = recalibrated_bam_basename,
          recalibration_report = GatherBqsrReports.output_bqsr_report,
          sequence_group_interval = subgroup,
          ref_dict = references.reference_fasta.ref_dict,
          ref_fasta = references.reference_fasta.ref_fasta,
          ref_fasta_index = references.reference_fasta.ref_fasta_index,
          bqsr_scatter = bqsr_divisor,
          compression_level = compression_level,
          preemptible_tries = papi_settings.agg_preemptible_tries,
          bin_base_qualities = bin_base_qualities,
          somatic = somatic
      }
    }

    # Merge the recalibrated BAM files resulting from by-interval recalibration
    call Processing.GatherSortedBamFiles as GatherBamFiles {
      input:
        input_bams = ApplyBQSR.recalibrated_bam,
        output_bam_basename = sample_and_bam.base_file_name,
        total_input_size = agg_bam_size,
        compression_level = compression_level,
        preemptible_tries = papi_settings.agg_preemptible_tries
    }
  }

  # Outputs that will be retained when execution is complete
  output {
    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File? selfSM = CheckContamination.selfSM
    Float? contamination = CheckContamination.contamination

    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File? output_bqsr_reports = GatherBqsrReports.output_bqsr_report

    File output_bam = SortSampleBam.output_bam
    File output_bam_index = select_first([SortSampleBam.output_bam_index])
  }
  meta {
    allowNestedInputs: true
  }
}

# create bam index
task IndexBam {
  input {
    File bam_input
  }

  # input file size
  Float input_size = size(bam_input, "GB")

  # output name for indexed bam
  String bam_index_output_name = bam_input + ".bai"

  command <<<
    samtools index -b ~{bam_input} ~{bam_index_output_name}
  >>>

  runtime {
    docker: "quay.io/broadinstitute/samtools:1.9"
    # if the input size is less than 1 GB adjust to min input size of 1 GB
    # disks should be set to 2 * input file size
    disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
    cpu: 1
    memory: "3.5 GB"
  }

  output {
    File bam_index_output = bam_index_output_name
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
