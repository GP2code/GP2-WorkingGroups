##  ExpansionHunter denovo (EHdn)
##
##  This WDL implements workflows for EHdn's individual profiling:
##  - case-control analysis (read EHdn docs on this:
##    https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/03_Case_control_quickstart.md);
##  - outlier analysis (read EHdn docs on this:
##    https://github.com/Illumina/ExpansionHunterDenovo/blob/master/documentation/04_Outlier_quickstart.md).

version 1.0

#import "Structs.wdl"
struct RuntimeAttr {
  Float? mem_gb
  Int? cpu_cores
  Int? disk_gb
  Int? boot_disk_gb
  Int? preemptible_tries
  Int? max_retries
}

struct FilenamePostfixes {
  String locus
  String motif
  String profile
  String merged_profile
  Int profile_len
}

workflow EHdnSTRprofile {

  input {
    File sample_bam_or_cram
    File sample_bam_or_cram_index
    File reference_fasta
    File reference_fasta_index
    Int min_anchor_mapq
    Int max_irr_mapq
    String ehdn_docker = "zihhuafang/expansionunterdenovo:v0.9.0"
    RuntimeAttr? runtime_attr_str_profile
    String runtime_zones
  }

  parameter_meta {
    sample_bam_or_cram: "bam or cram file to be used as input to EHdn."
    sample_bam_or_cram_index: "[Optional] index files for the sample bam/cram. Files should be in the same order as the samples. If not provided, it will be inferred from sample filenames."
    reference_fasta: "Sets the path to the reference."
    reference_fasta_index: "[Optional] Sets the path to the index of reference. If not provided, it will be inferred from the reference filename."
    min_anchor_mapq: "Sets anchor mapping quality (mapq) threshold."
    max_irr_mapq: "Set the mapping quality (mapq) threshold on reads to be considered when searching for in-repeat reads (IRR)."
    ehdn_docker: "Sets the docker image of EHdn."
    runtime_attr_str_profile: "[Optional] Override the default runtime attributes for STR profiling task."
  }

  # The values of the variables are based on
  # EHdn's current latest hard-coded postfixes.
  # Do not change them unless they are changed
  # in EHdn.
  FilenamePostfixes postfixes = object {
    locus: ".locus.tsv",
    motif: ".motif.tsv",
    profile: ".str_profile.json",
    merged_profile: ".multisample_profile.json",

    # This the length of `profile` postfix without
    # the filename extension. It is used to remove
    # postfix in order to extract sample name from
    # EHdn generated output.
    # e.g., extract `sample1` from `sample1.str_profile.json`.
    profile_len: 12
  }

  Boolean is_bam =
    basename(sample_bam_or_cram, ".bam") + ".bam" ==
    basename(sample_bam_or_cram)

  String filename =
    if is_bam then
      basename(sample_bam_or_cram, ".bam")
    else
      basename(sample_bam_or_cram, ".cram")

  call ComputeSTRProfile {
    input:
      filename = filename,
      bam_or_cram = sample_bam_or_cram,
      bam_or_cram_index = sample_bam_or_cram_index,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      min_anchor_mapq = min_anchor_mapq,
      max_irr_mapq = max_irr_mapq,
      ehdn_docker = ehdn_docker,
      runtime_attr_override = runtime_attr_str_profile,
      postfixes = postfixes,
      runtime_zones= runtime_zones
  }

  output {
    File sample_ehdn_locus = ComputeSTRProfile.locus
    File sample_ehdn_motif = ComputeSTRProfile.motif
    File sample_ehdn_str_profile = ComputeSTRProfile.str_profile
  }
}

task ComputeSTRProfile {
  input {
    String filename
    File bam_or_cram
    File bam_or_cram_index
    File reference_fasta
    File? reference_fasta_index
    Int min_anchor_mapq
    Int max_irr_mapq
    String ehdn_docker
    RuntimeAttr? runtime_attr_override
    FilenamePostfixes postfixes
    String runtime_zones
  }

  output {
    File locus = "${filename}${postfixes.locus}"
    File motif = "${filename}${postfixes.motif}"
    File str_profile = "${filename}${postfixes.profile}"
  }

  # Defining this varialbe is requires since it
  # will trigger localization of the file. The
  # localized file is used by the `profile` subcommand
  # of ExpansionHunterDenovo.
  File reference_fasta_index_ = select_first([
    reference_fasta_index,
    reference_fasta + ".fai"])

  command <<<
    set -euxo pipefail

    ExpansionHunterDenovo profile \
    --reads ~{bam_or_cram} \
    --reference ~{reference_fasta} \
    --output-prefix ~{filename} \
    --min-anchor-mapq ~{min_anchor_mapq} \
    --max-irr-mapq ~{max_irr_mapq}

  >>>

  RuntimeAttr runtime_attr_str_profile_default = object {
    cpu_cores: 1,
    mem_gb: 25,
    boot_disk_gb: 20,
    preemptible_tries: 2,
    max_retries: 1,
    disk_gb: 20 + ceil(size([
      bam_or_cram,
      bam_or_cram_index,
      reference_fasta,
      reference_fasta_index_], "GiB"))
  }
  RuntimeAttr runtime_attr = select_first([
    runtime_attr_override,
    runtime_attr_str_profile_default])

  runtime {
    docker: ehdn_docker
    cpu: runtime_attr.cpu_cores
    memory: runtime_attr.mem_gb + " GiB"
    disks: "local-disk " + runtime_attr.disk_gb + " HDD"
    bootDiskSizeGb: runtime_attr.boot_disk_gb
    preemptible: runtime_attr.preemptible_tries
    maxRetries: runtime_attr.max_retries
    zones: runtime_zones
  }
}