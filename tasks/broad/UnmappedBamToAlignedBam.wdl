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

# WORKFLOW DEFINITION
workflow UnmappedBamToAlignedBam {

  input {
    SampleAndBam sample_and_bam
    DNASeqSingleSampleReferences references
    PapiSettings papi_settings

    Boolean check_contamination = true
    Boolean check_fingerprints = true

    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu

    String cross_check_fingerprints_by
    File haplotype_database_file
    Float lod_threshold
    Boolean hard_clip_reads = false
    
    Boolean to_cram = false
    String? subset_region
  }

  Float cutoff_for_large_rg_in_gb = 7.5

  # Get the size of the standard reference files as well as the additional reference files needed for BWA

  File unmapped_bam = sample_and_bam.bam_or_cram
  File unmapped_index = select_first([sample_and_bam.bai_or_crai, sub(unmapped_bam, ".bam$", ".bai")])
  if (basename(unmapped_bam, ".bam") == basename(unmapped_bam)) {
    # Not a BAM suffix, so probably CRAM
    File unmapped_bam_index = select_first([sample_and_bam.bai_or_crai, sub(unmapped_bam, ".cram$", ".crai")])
  }

  # We don't partition input BAM for distributed alignment, to avoid extra copying
  # and extra deduplication+merge step with a large overhead.
  # However, if you want shared parallelilsm, consider adding 
  # bazam sharded parallelism like in this pipeline:
  # https://github.com/Oshlack/STRetch/blob/c5345e5dea4adfde790befb9903ec2d81ed5b2c1/pipelines/pipeline_stages.groovy#L101
  call Alignment.Bazam as Bazam {
    input:
      input_bam = unmapped_bam,
      input_bam_index = unmapped_index,
      sample_name = sample_and_bam.sample_name,
      output_bam_basename = sample_and_bam.base_file_name,
      reference_fasta = references.reference_fasta,
      preemptible_tries = papi_settings.preemptible_tries,
      duplicate_metrics_fname = sample_and_bam.base_file_name + ".duplicate_metrics",
      to_cram = to_cram,
      subset_region = subset_region,
  }

  File output_aligned_bam = Bazam.output_bam
  File output_aligned_bam_index = Bazam.output_bam_index

  Float mapped_bam_size = size(output_aligned_bam, "GiB")
  File mapped_bam = output_aligned_bam
  File mapped_bam_index = output_aligned_bam_index

  if (defined(haplotype_database_file) && check_fingerprints) {
    # Check identity of fingerprints across readgroups
    call QC.CrossCheckFingerprints as CrossCheckFingerprints {
      input:
        input_bams = [ mapped_bam ],
        input_bam_indexes = [ mapped_bam_index ],
        haplotype_database_file = haplotype_database_file,
        metrics_filename = sample_and_bam.base_file_name + ".crosscheck",
        total_input_size = mapped_bam_size,
        lod_threshold = lod_threshold,
        cross_check_by = cross_check_fingerprints_by,
        preemptible_tries = papi_settings.agg_preemptible_tries,
        ref_dict = references.reference_fasta.ref_dict,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index
    }
  }

  if (check_contamination) {
    # Estimate level of cross-sample contamination
    call Processing.CheckContamination as CheckContamination {
      input:
        input_bam = mapped_bam,
        input_bam_index = mapped_bam_index,
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

  # Outputs that will be retained when execution is complete
  output {
    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File? selfSM = CheckContamination.selfSM
    Float? contamination = CheckContamination.contamination

    File duplicate_metrics = Bazam.duplicate_metrics

    File output_bam = mapped_bam
    File output_index = mapped_bam_index
  }
  meta {
    allowNestedInputs: true
  }
}
