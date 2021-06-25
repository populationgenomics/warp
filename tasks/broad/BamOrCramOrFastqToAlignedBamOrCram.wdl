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
workflow BamOrCramOrFastqToAlignedBamOrCram {

  input {
    Input inp
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

  # We don't partition input BAM for distributed alignment, to avoid extra copying
  # and extra deduplication+merge step with a large overhead.
  # However, if you want shared parallelilsm, consider adding 
  # bazam sharded parallelism like in this pipeline:
  # https://github.com/Oshlack/STRetch/blob/c5345e5dea4adfde790befb9903ec2d81ed5b2c1/pipelines/pipeline_stages.groovy#L101

  if (inp.bam_or_cram_or_fastq1 == sub(inp.bam_or_cram_or_fastq1, ".cram$", "") + ".cram" ||
    inp.bam_or_cram_or_fastq1 == sub(inp.bam_or_cram_or_fastq1, ".bam$", "") + ".bam") {
  
    call Alignment.BwaFromBamOrCram {
      input:
        bam_or_cram = inp.bam_or_cram_or_fastq1,
        bai_or_crai = inp.bai_or_crai_or_fastq2,
        sample_name = inp.sample_name,
        output_bam_basename = inp.base_file_name,
        reference_fasta = references.reference_fasta,
        preemptible_tries = papi_settings.preemptible_tries,
        duplicate_metrics_fname = inp.base_file_name + ".duplicate_metrics",
        to_cram = to_cram,
        subset_region = subset_region,
    }
  }

  if (inp.bam_or_cram_or_fastq1 != sub(inp.bam_or_cram_or_fastq1, ".cram$", "") + ".cram" &&
    inp.bam_or_cram_or_fastq1 != sub(inp.bam_or_cram_or_fastq1, ".bam$", "") + ".bam") {
  
    call Alignment.BwaFromFastq {
      input:
        fastq1 = inp.bam_or_cram_or_fastq1,
        fastq2 = inp.bai_or_crai_or_fastq2,
        sample_name = inp.sample_name,
        output_bam_basename = inp.base_file_name,
        reference_fasta = references.reference_fasta,
        preemptible_tries = papi_settings.preemptible_tries,
        duplicate_metrics_fname = inp.base_file_name + ".duplicate_metrics",
        to_cram = to_cram
    }
  }
  
  File mapped_file = select_first([BwaFromBamOrCram.output_file, BwaFromFastq.output_file])
  File mapped_indx = select_first([BwaFromBamOrCram.output_indx, BwaFromFastq.output_indx])
  Float mapped_file_size = size(mapped_file, "GiB")

  if (defined(haplotype_database_file) && check_fingerprints) {
    # Check identity of fingerprints across readgroups
    call QC.CrossCheckFingerprints as CrossCheckFingerprints {
      input:
        input_bams = [ mapped_file ],
        input_bam_indexes = [ mapped_indx ],
        haplotype_database_file = haplotype_database_file,
        metrics_filename = inp.base_file_name + ".crosscheck",
        total_input_size = mapped_file_size,
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
        input_bam = mapped_file,
        input_bam_index = mapped_indx,
        contamination_sites_ud = contamination_sites_ud,
        contamination_sites_bed = contamination_sites_bed,
        contamination_sites_mu = contamination_sites_mu,
        ref_fasta = references.reference_fasta.ref_fasta,
        ref_fasta_index = references.reference_fasta.ref_fasta_index,
        output_prefix = inp.base_file_name,
        preemptible_tries = papi_settings.agg_preemptible_tries,
        contamination_underestimation_factor = 0.75
    }
  }

  # Outputs that will be retained when execution is complete
  output {
    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.cross_check_fingerprints_metrics

    File? selfSM = CheckContamination.selfSM
    Float? contamination = CheckContamination.contamination

    File duplicate_metrics = select_first([BwaFromBamOrCram.duplicate_metrics, BwaFromFastq.duplicate_metrics])

    File output_file = mapped_file
    File output_indx = mapped_indx
  }
  meta {
    allowNestedInputs: true
  }
}
