version 1.0

## Copyright Broad Institute, 2018
##
## This WDL defines tasks used for alignment of human whole-genome or exome sequencing data.
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

import "../../structs/dna_seq/DNASeqStructs.wdl"

task BwaFromFastq {
  input {
    File fastq1
    File fastq2
    String sample_name
    String output_bam_basename
    String duplicate_metrics_fname

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int preemptible_tries
    Boolean to_cram = false
  }
  
  String output_format = if to_cram then "cram" else "bam"
  
  Int bwa_cpu = 25
  Int bamsormadup_cpu = 6
  Int total_cpu = bwa_cpu + bamsormadup_cpu

  String rg_line = "@RG\\tID:~{sample_name}\\tSM:~{sample_name}"
  
  String output_file = if to_cram then "~{output_bam_basename}.cram" else "~{output_bam_basename}.bam"
  String output_indx = if to_cram then "~{output_bam_basename}.crai" else "~{output_bam_basename}.bai"
  
  # BWA command options:
  # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
  # -v3    minimum score to output [30]
  # -t16   threads
  # -Y     use soft clipping for supplementary alignments
  # -R     read group header line such as '@RG\tID:foo\tSM:bar'
  # -M     mark shorter split hits as secondary
  command <<<
    set -o pipefail
    set -ex

    (while true; do df -h; pwd; du -sh *; free -m; sleep 300; done) &
    
    bwa mem -K 100000000 -v3 -t~{bwa_cpu} -Y -R '~{rg_line}' \
      ~{reference_fasta.ref_fasta} ~{fastq1} ~{fastq2} \
      2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    bamsormadup inputformat=sam threads=~{bamsormadup_cpu} SO=coordinate \
      M=~{duplicate_metrics_fname} \
      outputformat=sam | \
    samtools view -T ~{reference_fasta.ref_fasta} \
      -O ~{output_format} \
      -o ~{output_file}
    
    samtools index -@~{total_cpu} ~{output_file} ~{output_indx}

    df -h; pwd; du -sh *
  >>>
  
  runtime {
    # docker: "australia-southeast1-docker.pkg.dev/fewgenomes/images/bazam:v2"
    # cromwell doesn't work with artifact registry:
    # java.lang.Exception: Registry australia-southeast1-docker.pkg.dev is not supported
    docker: "gcr.io/fewgenomes/bazam:v2"
    preemptible: preemptible_tries
    memory: "64 GiB"
    cpu: total_cpu
    disks: "local-disk " + 300 + " HDD"
  }
  output {
    File output_file = output_file
    File output_indx = output_indx
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{duplicate_metrics_fname}"
  }
}

task BwaFromBamOrCram {
  input {
    File bam_or_cram
    File bai_or_crai
    String sample_name
    String output_bam_basename
    String duplicate_metrics_fname

    # reference_fasta.ref_alt is the .alt file from bwa-kit
    # (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    ReferenceFasta reference_fasta

    Int preemptible_tries
    Boolean to_cram = false
    String? subset_region 
  }
  
  String output_format = if to_cram then "cram" else "bam"
  String bazam_regions = if defined(subset_region) then "--regions ~{subset_region} " else ""
  
  Int bwa_cpu = 20
  Int bazam_cpu = 5
  Int bamsormadup_cpu = 6
  Int total_cpu = bwa_cpu + bazam_cpu + bamsormadup_cpu

  String rg_line = "@RG\\tID:~{sample_name}\\tSM:~{sample_name}"
  
  String output_file = if to_cram then "~{output_bam_basename}.cram" else "~{output_bam_basename}.bam"
  String output_indx = if to_cram then "~{output_bam_basename}.crai" else "~{output_bam_basename}.bai"
  
  # BWA command options:
  # -K     process INT input bases in each batch regardless of nThreads (for reproducibility)
  # -p     smart pairing (ignoring in2.fq)
  # -v3    minimum score to output [30]
  # -t16   threads
  # -Y     use soft clipping for supplementary alignments
  # -R     read group header line such as '@RG\tID:foo\tSM:bar'
  # -M     mark shorter split hits as secondary
  command <<<
    set -o pipefail
    set -ex

    (while true; do df -h; pwd; du -sh *; free -m; sleep 300; done) &
    
    bazam -Xmx16g -Dsamjdk.reference_fasta=~{reference_fasta.ref_fasta} \
      ~{bazam_regions} -n~{bazam_cpu} -bam ~{bam_or_cram} | \
    bwa mem -K 100000000 -p -v3 -t~{bwa_cpu} -Y -R '~{rg_line}' \
      ~{reference_fasta.ref_fasta} /dev/stdin - \
      2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
    bamsormadup inputformat=sam threads=~{bamsormadup_cpu} SO=coordinate \
      M=~{duplicate_metrics_fname} \
      outputformat=sam | \
    samtools view -T ~{reference_fasta.ref_fasta} \
      -O ~{output_format} \
      -o ~{output_file}
    
    samtools index -@~{total_cpu} ~{output_file} ~{output_indx}

    df -h; pwd; du -sh *
  >>>
  
  runtime {
    # docker: "australia-southeast1-docker.pkg.dev/fewgenomes/images/bazam:v2"
    # cromwell doesn't work with artifact registry:
    # java.lang.Exception: Registry australia-southeast1-docker.pkg.dev is not supported
    docker: "gcr.io/fewgenomes/bazam:v2"
    preemptible: preemptible_tries
    memory: "64 GiB"
    cpu: total_cpu
    disks: "local-disk " + 300 + " HDD"
  }
  output {
    File output_file = output_file
    File output_indx = output_indx
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
    File duplicate_metrics = "~{duplicate_metrics_fname}"
  }
}

task SamSplitter {
  input {
    File input_bam
    Int n_reads
    Int preemptible_tries
    Int compression_level
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
  Float disk_multiplier = 2.5
  Int disk_size = ceil(disk_multiplier * unmapped_bam_size + 20)

  command {
    set -e
    mkdir output_dir

    total_reads=$(samtools view -c ~{input_bam})

    java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
      INPUT=~{input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=~{n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: preemptible_tries
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
}
