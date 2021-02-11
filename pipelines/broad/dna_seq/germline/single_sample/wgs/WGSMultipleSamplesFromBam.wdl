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
    Array[File] output_vcfs = select_all(WGSFromBam.output_vcf)
    Array[File] output_vcf_idx = select_all(WGSFromBam.output_vcf_index)
  }
  meta {
    allowNestedInputs: true
  }
}
