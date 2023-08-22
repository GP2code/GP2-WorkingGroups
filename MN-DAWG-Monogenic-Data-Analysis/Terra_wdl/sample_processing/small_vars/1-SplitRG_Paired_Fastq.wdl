version 1.0

# WORKFLOW DEFINITION

workflow SplitRGPairedFastqWF {
  input {
    Int threads
    String GP2_id 
    String fastq_R1 
    String fastq_R2 
    String platform_name
    String sequencing_center
    String run_date

    String docker_gdc = "quay.io/kmhernan/gdc-fastq-splitter:latest"
    String docker_fastp = "quay.io/biocontainers/fastp:0.20.1--h2e03b76_1"
    String runtime_zones 
  }
  # Run FASTP
  call fastp {
    input:
      output_prefix = GP2_id,
      threads = threads,
      fastq_R1 = fastq_R1,
      fastq_R2 = fastq_R2,
      docker_fastp = docker_fastp,
      runtime_zones =runtime_zones
  }
  # Split a fastq that has multiple readgroups
  call gdcfastqsplitter {
    input:
      output_prefix = GP2_id + "_",
      fastq_R1_qc = fastp.fastq_R1_qc,
      fastq_R2_qc = fastp.fastq_R2_qc,
      docker_gdc = docker_gdc,
      runtime_zones = runtime_zones
  }
    # Write data table file
  call datatablefile {
    input:
      output_prefix = GP2_id,
      fastq_r1 = gdcfastqsplitter.fastq_r1,
      fastq_r2 = gdcfastqsplitter.fastq_r2,
      platform_name = platform_name,
      sequencing_center = sequencing_center,
      run_date = run_date  
  }
  output {
    File fastq_report = fastp.fastq_report
    File readgroup_list = datatablefile.readgroup_list
  }
}

#__________________________________________________________________________________________________________

#TASK DEFINITIONS

# 1. run FASTP
task fastp {
  input {
    # Command parameters
    String output_prefix
    File fastq_R1
    File fastq_R2
    Int threads
    
    # Runtime parameters
    Int additional_disk_space_gb = 10
    Int machine_mem_gb = 10
    Int preemptible_attempts = 2
    String docker_fastp
    String runtime_zones 
  }
    Int command_mem_gb = machine_mem_gb - 1
    Int disk_space_gb = ceil((size(fastq_R1, "GB") + size(fastq_R2, "GB")) * 2 ) + additional_disk_space_gb
  command {
    fastp --thread ~{threads} -i ~{fastq_R1} -o ~{output_prefix}_R1_qc.fastq.gz -I ~{fastq_R2} -O ~{output_prefix}_R2_qc.fastq.gz \
                -j ~{output_prefix}.JSON -q 15 \
                -u 40 \
                -g >/dev/null
  }
  runtime {
    docker: docker_fastp
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    cpu: threads
    preemptible: preemptible_attempts
    zones: runtime_zones 
  }
  output {
    File fastq_R1_qc = "~{output_prefix}_R1_qc.fastq.gz"
    File fastq_R2_qc = "~{output_prefix}_R2_qc.fastq.gz"
    File fastq_report = "~{output_prefix}.JSON"
  }
}  
  
  
# 2. Split a fastq that has multiple readgroups
task gdcfastqsplitter {
  input {
    # Command parameters
    String output_prefix
    File fastq_R1_qc
    File fastq_R2_qc

    # Runtime parameters
    Int additional_disk_space_gb = 10
    Int machine_mem_gb = 10
    Int preemptible_attempts = 2
    String docker_gdc
    String runtime_zones 
  }
    Int command_mem_gb = machine_mem_gb - 1
    Int disk_space_gb = ceil((size(fastq_R1_qc, "GB") + size(fastq_R2_qc, "GB")) * 2 ) + additional_disk_space_gb
  command {
    gdc-fastq-splitter --output-prefix ~{output_prefix} ~{fastq_R1_qc} ~{fastq_R2_qc}
  }
  runtime {
    docker: docker_gdc
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    zones: runtime_zones 
  }
  output {
    Array[File] fastq_r1 = glob("*_R1.fq.gz")
    Array[File] fastq_r2 = glob("*_R2.fq.gz")   
  }
}


# 3. Write names of splitted fastq files into a txt file
task datatablefile {
  input {
    # Command parameters
    String output_prefix
    Array[String] fastq_r1
    Array[String] fastq_r2
    String library_name
    String platform_name
    String sequencing_center
    String run_date
    
    # Output tsv without header
    # Header is sample_name readgroup fastq1 fastq2 platform_unit platform_name library_name run_date sequencing_center
  }
  command <<<
    set -oe pipefail
    
    python << CODE
    
    import os
    filepaths1 = ['~{sep="','" fastq_r1}']
    filepaths2 = ['~{sep="','" fastq_r2}']
    platform_name  = '~{platform_name}'
    sequencing_center = '~{sequencing_center}'
    run_date = '~{run_date}'

    with open("~{output_prefix}.txt", "w") as fi:
        for i in range(len(filepaths1)):
            sample_name, flowcell, lane, suffix = os.path.split(filepaths1[i])[1].rsplit("_", 3)
            fi.write(sample_name + "\t" + flowcell + "_" + lane + "\t" + filepaths1[i] + "\t" + filepaths2[i] + "\t" + flowcell + "_" + lane + "\t" + platform_name + "\t" + sample_name  + "\t" + run_date + "\t" + sequencing_center + "\n" )
    CODE
    
    >>>
  output {
    File readgroup_list = "~{output_prefix}.txt"
  }
  runtime {
    docker: "python:latest"
    preemptible: 3
  }
}