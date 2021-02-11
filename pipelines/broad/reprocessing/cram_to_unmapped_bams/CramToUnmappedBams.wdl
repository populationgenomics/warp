version 1.0

# Exactly one of input_cram and input_bam should be supplied to this workflow. If an input_cram is supplied, a ref_fasta
# and ref_fasta_index must also be supplied. The ref_fasta and ref_fasta_index are used to generate a bam, so if an
# input_cram is not supplied, the input_bam is used instead and the ref_fasta and ref_fasta_index are not needed.

# If the output_map file is provided, it is expected to be a tab-separated file containing a list of all the read group ids
# found in the input_cram / input_bam and the desired name of the unmapped bams generated for each.
# If the file is not provided, the output names of the unmapped bams will be the read_group_id<unmapped_bam_suffix>
workflow BamToFastq {

  String pipeline_version = "1.0.0"

  input {
    File? input_cram
    File? input_bam
    File? ref_fasta
    File? ref_fasta_index
    Int additional_disk = 20
  }

  File input_file = input_bam
  Float input_size = size(input_file, "GiB")

  String output_fq_basename = basename(input_file, ".bam")

  output {
    File unmapped_bam = SortSam.output_bam
  }
  meta {
    allowNestedInputs: true
  }
}

task SubsetBam {
  input {
    File bam
    String output_base_name
  }

  Int disk_size = ceil(size(bam, "GiB")) + 10

  command {
    samtools index ~{bam}
    samtools view ~{bam} 21 -Obam -o ~{output_base_name}.bam
    samtools index ~{output_base_name}.bam
    mv ~{output_base_name}.bam.bai ~{output_base_name}.bai
  }

  output {
    File output_bam = "~{output_base_name}.bam"
    File output_bam_index = "~{output_base_name}.bai"
  }

  runtime {
    docker: "biocontainers/samtools:1.3.1"
    memory: "4 GiB"
    disks: "local-disk ~{disk_size} HDD"
    preemptible: 3
    cpu: 1
  }
}

task RevertSam {
  input {
    File input_bam
    String output_bam_filename
    Int disk_size
    Int memory = 3
  }

  command <<<
    java -Xms3500m -jar /usr/picard/picard.jar \
    RevertSam \
    --INPUT ~{input_bam} \
    --OUTPUT ~{output_bam_filename} \
    --VALIDATION_STRINGENCY LENIENT \
    --ATTRIBUTE_TO_CLEAR FT \
    --ATTRIBUTE_TO_CLEAR CO \
    --ATTRIBUTE_TO_CLEAR PA \
    --ATTRIBUTE_TO_CLEAR OA \
    --ATTRIBUTE_TO_CLEAR XA \
    --SORT_ORDER coordinate

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory} GiB"
    preemptible: 3
  }

  output {
    File output_bam = output_bam_filename
  }
}

# This task is slower than converting straight from cram to bam (note we stream through sam format in between cram and bam)
# This is currently necessary due to a difference in the way the NM tag is calculated in order to produce a valid bam.
task CramToBam {
  input {
    File ref_fasta
    File ref_fasta_index
    File cram_file
    String output_basename
    Int disk_size
    Int memory = 7
  }

  command <<<

    set -e
    set -o pipefail

    samtools view -h -T ~{ref_fasta} ~{cram_file} |
    samtools view -b -o ~{output_basename}.bam -
    samtools index -b ~{output_basename}.bam
    mv ~{output_basename}.bam.bai ~{output_basename}.bai

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    cpu: 3
    memory: "~{memory} GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
  }

  output {
    File output_bam = "~{output_basename}.bam"
    File output_bam_index = "~{output_basename}.bai"
  }
}

task GenerateOutputMap {
  input {
    File input_bam
    String unmapped_bam_suffix
    Int disk_size
    Int memory = 3
  }

  command <<<

    set -e

    samtools view -H ~{input_bam} | grep '^@RG' | cut -f2 | sed s/ID:// > readgroups.txt

    echo -e "#READ_GROUP_ID\tOUTPUT" > output_map.tsv

    for rg in `cat readgroups.txt`; do
      echo -e "$rg\t$rg~{unmapped_bam_suffix}" >> output_map.tsv
    done

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory} GiB"
    preemptible: 3
  }

  output {
    File output_map = "output_map.tsv"
  }
}

task SplitUpOutputMapFile {
  input {
    File read_group_map_file
    Int disk_size = 10
    Int memory = 3
  }

  command <<<
    mkdir rgtemp
    cd rgtemp

    # splits list of mappings into single files.  One line each.
    grep -v '^#' ~{read_group_map_file} | split -l 1 - rg_to_ubam_
  >>>

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory} GiB"
  }

  output {
    Array[File] rg_to_ubam_file = glob("rgtemp/rg_to_ubam_*")
  }
}

task SplitOutUbamByReadGroup {

  input {
    File input_bam
    File rg_to_ubam_file
    Int disk_size
    Int memory = 30
  }

  Array[Array[String]] tmp = read_tsv(rg_to_ubam_file)

  command <<<
    echo "Read Group ~{tmp[0][0]} from ~{input_bam} is being written to ~{tmp[0][1]}"
    samtools view -b -h -r ~{tmp[0][0]} -o ~{tmp[0][1]} ~{input_bam}
  >>>

  output {
    File output_bam = tmp[0][1]
  }

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    cpu: 2
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory} GiB"
    preemptible: 3
  }
}

task ValidateSamFile {
  input {
    File input_bam
    String report_filename
    Int disk_size
    Int memory = 3
  }

  command <<<

    java -Xms3500m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      --INPUT ~{input_bam} \
      --OUTPUT ~{report_filename} \
      --MODE VERBOSE \
      --IS_BISULFITE_SEQUENCED false

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory} GiB"
    preemptible: 3
  }

  output {
    File report = "~{report_filename}"
  }
}

task SortSam {
  input {
    File input_bam
    String output_bam_filename
    Int memory = 7
    Float sort_sam_disk_multiplier = 6
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20

  command <<<
    java -Xms7g -jar /usr/picard/picard.jar \
    SortSam \
    --INPUT ~{input_bam} \
    --OUTPUT ~{output_bam_filename} \
    --SORT_ORDER queryname \
    --MAX_RECORDS_IN_RAM 300000

  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    disks: "local-disk " + disk_size + " HDD"
    memory: "~{memory} GiB"
    preemptible: 3
  }

  output {
    File output_bam = output_bam_filename
  }
}
