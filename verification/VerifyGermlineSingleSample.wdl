version 1.0

#import "../verification/VerifyMetrics.wdl" as MetricsVerification
#import "../verification/VerifyTasks.wdl" as Tasks

workflow VerifyGermlineSingleSample {

  input {
    Array[File] truth_metrics
    Array[File] test_metrics

    File truth_cram
    File truth_crai
    File test_cram
    File test_crai

    File truth_gvcf
    File test_gvcf
  }

#  call MetricsVerification.VerifyMetrics as CompareMetrics {
#    input:
#      test_metrics = test_metrics,
#      truth_metrics = truth_metrics
#  }

  call CompareGvcfs {
    input:
      test_gvcf = test_gvcf,
      truth_gvcf = truth_gvcf
  }

  call SummarizeResults {
    input:
      compare_gvcfs_exit_code = CompareGvcfs.exit_code,
      compare_gvcfs_results_file = CompareGvcfs.error_report_file
  }

#  call Tasks.CompareCrais {
#    input:
#      test_crai = test_crai,
#      truth_crai = truth_crai
#  }
#
#  call Tasks.CompareCrams {
#    input:
#      test_cram = test_cram,
#      test_crai = test_crai,
#      truth_cram = truth_cram,
#      truth_crai = truth_crai
#  }

  output {
    File report_file = SummarizeResults.error_report_file
#    Array[File] metric_comparison_report_files = CompareMetrics.metric_comparison_report_files
  }
  meta {
    allowNestedInputs: true
  }
}


task SummarizeResults {
  input {
    Int compare_gvcfs_exit_code
    File compare_gvcfs_results_file
  }

  command {
    cat ~{compare_gvcfs_results_file}
  }
  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    File error_report_file = stdout()
  }
}

task CompareGvcfs {
  input {
    File test_gvcf
    File truth_gvcf
  }

  Int ret_val = 0

  command {
    echo "Results of CompareGvcfs:"
    echo -e "Test:\t~{test_gvcf}"
    echo -e "Truth:\t~{truth_gvcf}"
    DIFF_LINES=$(diff <(gunzip -c -f ~{test_gvcf} | grep -v '^##') <(gunzip -c -f ~{truth_gvcf} | grep -v '^##') | grep -e "^<" | wc -l)
    if [ $DIFF_LINES -lt 10 ]; then
      echo -e "Pass\tGVCs differ by $DIFF_LINES lines"
    else
      ~{ret_val}=1
      echo "Fail\tGVCFs differ by $DIFF_LINES lines"
      DIFF_LINES=$(diff <(gunzip -c -f ~{test_gvcf} | grep -v '^##' | cut -f 1-5,7-) <(gunzip -c -f ~{truth_gvcf} | grep -v '^##' | cut -f 1-5,7-) | grep -e "^<" | wc -l)
      if [ $DIFF_LINES -eq 0 ]; then
        echo -e "\tNote that differences are ONLY found in the 'quality' column"
      fi
    fi
    exit 0
  }

  runtime {
    docker: "gcr.io/gcp-runtimes/ubuntu_16_0_4:latest"
    disks: "local-disk 50 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Int exit_code = ret_val
    File error_report_file = stdout()
  }
}
