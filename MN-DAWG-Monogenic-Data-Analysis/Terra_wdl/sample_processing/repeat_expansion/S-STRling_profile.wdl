##  STRling
##
##  This WDL implements workflow for STRling's individual profiling:

version 1.0

struct RuntimeAttr {
  Float? mem_gb
  Int? cpu_cores
  Int? disk_gb
  Int? boot_disk_gb
  Int? preemptible_tries
  Int? max_retries
}

workflow strling_profile {
  input {
    File sample_bam_or_cram
    File sample_bam_or_cram_index
    File reference_fasta
    File reference_fasta_index
    String? strling_docker_override
    String strling_docker = select_first([strling_docker_override, "quay.io/biocontainers/strling:0.5.1--h14cfee4_0"])
    RuntimeAttr? runtime_attr_str_profile   
    String runtime_zones
  }
  
  parameter_meta {
    sample_bam_or_cram: "bam or cram file to be used as input to STRling."
    sample_bam_or_cram_index: "[Optional] index files for the sample bam/cram. Files should be in the same order as the samples. If not provided, it will be inferred from sample filenames."
    reference_fasta: "Sets the path to the reference."
    reference_fasta_index: "[Optional] Sets the path to the index of reference. If not provided, it will be inferred from the reference filename."
    strling_docker: "Sets the docker image of STRling."
    runtime_attr_str_profile: "[Optional] Override the default runtime attributes for STR profiling task."
  }

  Boolean is_bam =
    basename(sample_bam_or_cram, ".bam") + ".bam" ==
    basename(sample_bam_or_cram)

  String filename =
    if is_bam then
      basename(sample_bam_or_cram, ".bam")
    else
      basename(sample_bam_or_cram, ".cram")

  call str_extract {
    input:
      ref_fasta = reference_fasta,
      sampleID = filename,
      cram = sample_bam_or_cram,
      crai = sample_bam_or_cram_index,
      strling_docker = strling_docker,
      runtime_attr_override = runtime_attr_str_profile,
      runtime_zones = runtime_zones
  }
  
  output {
    File sample_STRling_bin = str_extract.bin
  }
}

task str_extract {
  input {
    File ref_fasta
    File cram
    File crai
    String sampleID
    RuntimeAttr? runtime_attr_override
    String runtime_zones
    String strling_docker
  }

  output {
    File bin = "~{sampleID}.bin"
  }

  command <<<
    strling extract \
      -f ~{ref_fasta} \
      ~{cram} \
      ~{sampleID}.bin
  >>>

  RuntimeAttr runtime_attr_str_profile_default = object {    
    cpu_cores: 1,
    mem_gb: 8,
    boot_disk_gb: 20,
    preemptible_tries: 2,
    max_retries: 1,
    disk_gb: 20 + ceil(size([
      cram,
      crai,
      ref_fasta], "GiB"))
  }
  
  RuntimeAttr runtime_attr = select_first([
    runtime_attr_override,
    runtime_attr_str_profile_default])

  runtime {
    docker: strling_docker
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: runtime_attr.boot_disk_gb
    preemptible: runtime_attr.preemptible_tries
    maxRetries: runtime_attr.max_retries    
    zones: runtime_zones
  }
}