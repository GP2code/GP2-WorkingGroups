version 1.0

# WORKFLOW DEFINITION

workflow fastq_to_cram {
  input {
    # fastq processing
    String fastq_R1 
    String fastq_R2 
    String platform_name
    String sequencing_center
    String run_date

    # sample-level info
    String sample_name 
    String unmapped_bam_suffix
    Int read_length

    #alignment QC files
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File? fingerprint_genotypes_file
    File? haplotype_database_file
    File wgs_evaluation_interval_list
    File wgs_coverage_interval_list
    File wgs_calling_interval_list

    # reference genome files
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File ref_alt
    File ref_bwt
    File ref_sa
    File ref_amb
    File ref_ann
    File ref_pac
    File somalier_sites

    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices

    # Optional input to increase all disk sizes in case of outlier sample with strange size behavior
    Int? increase_disk_size

    # MarkDuplicates and SortSam currently take too long for preemptibles if the input data is too large
    Float gb_size_cutoff_for_preemptibles = 110.0
    Boolean data_too_large_for_preemptibles = SumFloats.total_size > gb_size_cutoff_for_preemptibles

    # Some tasks need wiggle room, and we also need to add a small amount of disk to prevent getting a
    # Cromwell error from asking for 0 disk when the input is less than 1GB
    Int additional_disk = select_first([increase_disk_size, 20])
    
    # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation
    # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
    Float max_duplication_in_reasonable_sample = 0.30
    Float max_chimerism_in_reasonable_sample = 0.15
    
    String bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"
    Int compression_level = 5
    
    #dockers 
    Int? fastp_threads_override
    Int fastp_threads = select_first([fastp_threads_override, 4])
    String docker_gdc = "quay.io/kmhernan/gdc-fastq-splitter:latest"
    String docker_fastp = "quay.io/biocontainers/fastp:0.23.4--hadf994f_2"
    String? gatk_docker_override
    String gatk_docker = select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.2.6.0"])
    String? gatk_path_override
    String gatk_path = select_first([gatk_path_override, "gatk"])
    String gotc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String gotc_path = "/usr/gitc/"
    String python_docker = "us.gcr.io/broad-gotc-prod/python:2.7" 
    String somalier_docker = "brentp/somalier:v0.2.19"

    Int preemptible_tries
    Int agg_preemptible_tries
    Int max_retries = 1
    String? runtime_zones_override
    String runtime_zones = select_first([runtime_zones_override,"us-central1-a us-central1-b us-central1-c us-central1-f"])
  }


  # Run FASTP
  call fastp {
    input:
      output_prefix = sample_name,
      threads = fastp_threads,
      fastq_R1 = fastq_R1,
      fastq_R2 = fastq_R2,
      docker_fastp = docker_fastp,
      runtime_zones =runtime_zones
  }

  # Split a fastq that has multiple readgroups
  call gdcfastqsplitter {
    input:
      output_prefix = sample_name + "_",
      fastq_R1_qc = fastp.fastq_R1_qc,
      fastq_R2_qc = fastp.fastq_R2_qc,
      docker_gdc = docker_gdc,
      runtime_zones = runtime_zones
  }

  # Write data table file
  call datatablefile {
    input:
      output_prefix = sample_name,
      fastq_r1 = gdcfastqsplitter.fastq_r1,
      fastq_r2 = gdcfastqsplitter.fastq_r2,
      platform_name = platform_name,
      sequencing_center = sequencing_center,
      run_date = run_date  
  }

  Array[Array[String]] readgroup_array = read_tsv(datatablefile.readgroup_list) 

  # Convert multiple pairs of input fastqs in parallel
  scatter (i in range(length(readgroup_array))) {

    # Convert pair of FASTQs to uBAM
    call PairedFastQsToUnmappedBAM {
      input:
        sample_name = readgroup_array[i][0],
        readgroup_name = readgroup_array[i][1],
        fastq_1 = readgroup_array[i][2],
        fastq_2 = readgroup_array[i][3],
        platform_unit = readgroup_array[i][4],
        platform_name = readgroup_array[i][5],
        library_name = readgroup_array[i][6],
        run_date = readgroup_array[i][7],
        sequencing_center = readgroup_array[i][8],
        gatk_path = gatk_path,
        gatk_docker = gatk_docker,
        preemptible_attempts = preemptible_tries,
        runtime_zones = runtime_zones
    }
  }
  
  Array[File] flowcell_unmapped_bams = PairedFastQsToUnmappedBAM.output_unmapped_bam

  # Get the version of BWA to include in the PG record in the header of the BAM produced
  # by MergeBamAlignment.
  call GetBwaVersion {
    input: 
      gotc_docker = gotc_docker,
      bwa_path = gotc_path,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  # Align flowcell-level unmapped input bams in parallel
  scatter (unmapped_bam in flowcell_unmapped_bams) {

    Float unmapped_bam_size = size(unmapped_bam, "GB")

    String sub_strip_path = "gs://.*/"
    String sub_strip_unmapped = unmapped_bam_suffix + "$"
    String sub_sub = sub(sub(unmapped_bam, sub_strip_path, ""), sub_strip_unmapped, "")

    # QC the unmapped BAM
    call CollectQualityYieldMetrics {
      input:
        input_bam = unmapped_bam,
        metrics_filename = sub_sub + ".unmapped.quality_yield_metrics",
        preemptible_tries = preemptible_tries,
        runtime_zones = runtime_zones
    }

    # Map reads to reference
    call SamToFastqAndBwaMemAndMba {
      input:
        input_bam = unmapped_bam,
        bwa_commandline = bwa_commandline,
        output_bam_basename = sub_sub + ".aligned.unsorted",
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        ref_alt = ref_alt,
        ref_bwt = ref_bwt,
        ref_amb = ref_amb,
        ref_ann = ref_ann,
        ref_pac = ref_pac,
        ref_sa = ref_sa,
        bwa_version = GetBwaVersion.version,
        compression_level = compression_level,
        preemptible_tries = preemptible_tries,
        gotc_docker = gotc_docker,
        bwa_path = gotc_path,
        gotc_path = gotc_path,
        max_retries = max_retries,
        runtime_zones = runtime_zones
    }

    Float mapped_bam_size = size(SamToFastqAndBwaMemAndMba.output_bam, "GB")

    # QC the aligned but unsorted readgroup BAM
    # no reference as the input here is unsorted, providing a reference would cause an error
    call CollectUnsortedReadgroupBamQualityMetrics {
      input:
        input_bam = SamToFastqAndBwaMemAndMba.output_bam,
        output_bam_prefix = sub_sub + ".readgroup",
        preemptible_tries = preemptible_tries,
        runtime_zones = runtime_zones
    }
  }

  # Sum the read group bam sizes to approximate the aggregated bam size
  call SumFloats {
    input:
      sizes = mapped_bam_size,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones
  }

  # Aggregate aligned+merged flowcell BAM files and mark duplicates
  # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
  # to avoid having to spend time just merging BAM files.
  call MarkDuplicates {
    input:
      input_bams = SamToFastqAndBwaMemAndMba.output_bam,
      output_bam_basename = sample_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = sample_name + ".duplicate_metrics",
      # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs
      # and the merged output.
      compression_level = compression_level,
      preemptible_tries = if data_too_large_for_preemptibles then 0 else agg_preemptible_tries,
      gatk_docker = gatk_docker,
      gatk_path = gatk_path,
      total_input_size = SumFloats.total_size,
      runtime_zones = runtime_zones
  }

  Float agg_bam_size = size(MarkDuplicates.output_bam, "GB")

  # Sort aggregated+deduped BAM file and fix tags
  call SortSam as SortSampleBam {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = sample_name + ".aligned.duplicate_marked.sorted",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      # This task spills to disk so we need space for the input bam, the output bam, and any spillage.
      compression_level = compression_level,
      preemptible_tries = agg_preemptible_tries,
      gatk_docker = gatk_docker,
      gatk_path = gatk_path,
      runtime_zones = runtime_zones
  }

  if (defined(haplotype_database_file)) {
    # Check identity of fingerprints across readgroups
    call CrossCheckFingerprints {
      input:
        input_bams = [SortSampleBam.output_bam],
        input_bam_indexes = [SortSampleBam.output_bam_index],
        haplotype_database_file = haplotype_database_file,
        metrics_filename = sample_name + ".crosscheck",
        preemptible_tries = agg_preemptible_tries,
        total_input_size = SumFloats.total_size,
        runtime_zones = runtime_zones
    }
  }

  # Estimate level of cross-sample contamination
  call CheckContamination {
    input:
      input_bam = SortSampleBam.output_bam, 
      input_bam_index = SortSampleBam.output_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_prefix = sample_name + ".preBqsr",
      preemptible_tries = agg_preemptible_tries,
      contamination_underestimation_factor = 0.75,
      runtime_zones = runtime_zones
  }

  # QC the final BAM (consolidated after scattered BQSR)
  call CollectReadgroupBamQualityMetrics {
    input:
      input_bam = SortSampleBam.output_bam, 
      input_bam_index = SortSampleBam.output_bam_index,
      output_bam_prefix = sample_name + ".readgroup",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = agg_preemptible_tries,
      runtime_zones = runtime_zones
  }

  # QC the final BAM some more (no such thing as too much QC)
  call CollectAggregationMetrics {
    input:
      input_bam = SortSampleBam.output_bam, 
      input_bam_index = SortSampleBam.output_bam_index,
      output_bam_prefix = sample_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = agg_preemptible_tries,
      runtime_zones = runtime_zones
  }

  if (defined(haplotype_database_file) && defined(fingerprint_genotypes_file)) {
    # Check the sample BAM fingerprint against the sample array
    call CheckFingerprint {
      input:
        input_bam = SortSampleBam.output_bam, 
        input_bam_index = SortSampleBam.output_bam_index,
        haplotype_database_file = haplotype_database_file,
        genotypes = fingerprint_genotypes_file,
        output_basename = sample_name,
        sample = sample_name,
        preemptible_tries = agg_preemptible_tries,
        runtime_zones = runtime_zones
    }
  }

  # QC the sample WGS metrics (stringent thresholds)
  call CollectWgsMetrics {
    input:
      input_bam = SortSampleBam.output_bam, 
      input_bam_index = SortSampleBam.output_bam_index,
      metrics_filename = sample_name + ".wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      read_length = read_length,
      preemptible_tries = agg_preemptible_tries,
      runtime_zones = runtime_zones
  }

  # QC the sample raw WGS metrics (common thresholds)
  call CollectRawWgsMetrics {
    input:
      input_bam = SortSampleBam.output_bam, 
      input_bam_index = SortSampleBam.output_bam_index,
      metrics_filename = sample_name + ".raw_wgs_metrics",
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      wgs_coverage_interval_list = wgs_coverage_interval_list,
      read_length = read_length,
      preemptible_tries = agg_preemptible_tries,
      runtime_zones = runtime_zones
  }

  String sort_bam_basename = sample_name + ".aligned.duplicate_marked.sorted"

  # Generate a checksum per readgroup in the final BAM
  call CalculateReadGroupChecksum {
    input:
      input_bam = SortSampleBam.output_bam, 
      input_bam_index = SortSampleBam.output_bam_index,
      read_group_md5_filename = sort_bam_basename + ".bam.read_group_md5",
      preemptible_tries = agg_preemptible_tries,
      runtime_zones = runtime_zones
  }

  # Convert the final merged recalibrated BAM file to CRAM format
  call ConvertToCram {
    input:
      input_bam = SortSampleBam.output_bam, 
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_basename = sample_name,
      preemptible_tries = agg_preemptible_tries,
      runtime_zones = runtime_zones
  }

  Float cram_size = size(ConvertToCram.output_cram, "GB")

  # Check whether the data has massively high duplication or chimerism rates
  call CheckPreValidation {
    input:
      duplication_metrics = MarkDuplicates.duplicate_metrics,
      chimerism_metrics = CollectAggregationMetrics.alignment_summary_metrics,
      max_duplication_in_reasonable_sample = max_duplication_in_reasonable_sample,
      max_chimerism_in_reasonable_sample = max_chimerism_in_reasonable_sample,
      preemptible_tries = agg_preemptible_tries,
      runtime_zones = runtime_zones
 }

  # Validate the CRAM file
  call ValidateSamFile as ValidateCram {
    input:
      input_bam = ConvertToCram.output_cram,
      input_bam_index = ConvertToCram.output_cram_index,
      report_filename = sample_name + ".cram.validation_report",
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      ignore = ["MISSING_TAG_NM"],
      max_output = 1000000000,
      is_outlier_data = CheckPreValidation.is_outlier_data,
      preemptible_tries = agg_preemptible_tries,
      max_retries = max_retries,
      runtime_zones = runtime_zones
  }

  # Quick QC by somalier
  call somalier {
    input:
      input_bam = SortSampleBam.output_bam, 
      input_bam_index = SortSampleBam.output_bam_index,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      somalier_sites = somalier_sites,
      sample_name = sample_name,
      preemptible_tries = preemptible_tries,
      runtime_zones = runtime_zones,
      somalier_docker = somalier_docker
  }


  output {
    File fastq_report = fastp.fastq_report

    Array[File] quality_yield_metrics = CollectQualityYieldMetrics.metrics

    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf
    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics
    Array[File] unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf
    Array[File] unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics
    Array[File] unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf
    Array[File] unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics
    Array[File] unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf
    Array[File] unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics

    File read_group_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics
    File read_group_gc_bias_detail_metrics = CollectReadgroupBamQualityMetrics.gc_bias_detail_metrics
    File read_group_gc_bias_pdf = CollectReadgroupBamQualityMetrics.gc_bias_pdf
    File read_group_gc_bias_summary_metrics = CollectReadgroupBamQualityMetrics.gc_bias_summary_metrics

    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.metrics

    String outlier = CheckPreValidation.is_outlier_data 

    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination
    
    File calculate_read_group_checksum_md5 = CalculateReadGroupChecksum.md5_file

    File agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics
    File agg_bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics
    File agg_bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics
    File agg_gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
    File agg_gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
    File agg_gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics
    File agg_insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf
    File agg_insert_size_metrics = CollectAggregationMetrics.insert_size_metrics
    File agg_pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics
    File agg_pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics
    File agg_quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
    File agg_quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics

    File? fingerprint_summary_metrics = CheckFingerprint.summary_metrics
    File? fingerprint_detail_metrics = CheckFingerprint.detail_metrics

    File wgs_metrics = CollectWgsMetrics.metrics
    File raw_wgs_metrics = CollectRawWgsMetrics.metrics
    File duplicate_metrics = MarkDuplicates.duplicate_metrics
    File sorted_bam = SortSampleBam.output_bam
    File sorted_bam_index = SortSampleBam.output_bam_index

    File output_cram = ConvertToCram.output_cram
    File output_cram_index = ConvertToCram.output_cram_index
    File output_cram_md5 = ConvertToCram.output_cram_md5
    File validate_cram_file_report = ValidateCram.report
    File somalier_output = somalier.somalier_output
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
    
    # Runtime parameters (n1-standard-4)
    Int additional_disk_space_gb = 10
    Int machine_mem_gb = 15
    Int preemptible_attempts = 3
    String docker_fastp
    String runtime_zones 
  }
    Int disk_space_gb = ceil((size(fastq_R1, "GiB") + size(fastq_R2, "GiB")) * 2 ) + additional_disk_space_gb
  command <<<
    fastp --thread ~{threads} -i ~{fastq_R1} -o ~{output_prefix}_R1_qc.fastq.gz -I ~{fastq_R2} -O ~{output_prefix}_R2_qc.fastq.gz \
                -j ~{output_prefix}.JSON -q 15 \
                -u 40 \
                -g >/dev/null
  >>>
  runtime {
    docker: docker_fastp
    memory: machine_mem_gb + " GiB"
    disks: "local-disk " + disk_space_gb + " SSD"
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
    Int preemptible_attempts = 3
    String docker_gdc
    String runtime_zones 
  }
    Int disk_space_gb = ceil((size(fastq_R1_qc, "GB") + size(fastq_R2_qc, "GB")) * 2 ) + additional_disk_space_gb

  command <<<
    gdc-fastq-splitter --output-prefix ~{output_prefix} ~{fastq_R1_qc} ~{fastq_R2_qc}
  >>>
  runtime {
    docker: docker_gdc
    cpu: 1
    disks: "local-disk " + disk_space_gb + " SSD"
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
    cpu: 1
    preemptible: 3
  }
}

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  input {
    # Command parameters
    String sample_name
    File fastq_1
    File fastq_2
    String readgroup_name
    String library_name
    String platform_unit
    String run_date
    String platform_name
    String sequencing_center

    # Runtime parameters
    Int? additional_disk_gb_override
    Int machine_mem_gb = 7
    Int preemptible_attempts
    String gatk_docker
    String gatk_path
    String runtime_zones
  }
    Int command_mem_gb = machine_mem_gb - 1
    Int additional_disk_space_gb = select_first([additional_disk_gb_override, 25])
    Int disk_space_gb = ceil((size(fastq_1, "GiB") + size(fastq_2, "GiB")) * 4 ) + additional_disk_space_gb

  command <<<
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    FastqToSam \
    --FASTQ ~{fastq_1} \
    --FASTQ2 ~{fastq_2} \
    --OUTPUT ~{readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ~{readgroup_name} \
    --SAMPLE_NAME ~{sample_name} \
    --LIBRARY_NAME ~{library_name} \
    --PLATFORM_UNIT ~{platform_unit} \
    --RUN_DATE ~{run_date} \
    --PLATFORM ~{platform_name} \
    --SEQUENCING_CENTER ~{sequencing_center} 
  >>>
  runtime {
    docker: gatk_docker
    memory: machine_mem_gb + " GiB"
    cpu: "1"
    disks: "local-disk " + disk_space_gb + " SSD"
    preemptible: preemptible_attempts
    zones: runtime_zones
  }
  output {
    File output_unmapped_bam = "~{readgroup_name}.unmapped.bam"
  }
}

# Collect sequencing yield quality metrics
task CollectQualityYieldMetrics {
  input {
    File input_bam
    String metrics_filename
    Int preemptible_tries
    String runtime_zones
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 20

  command <<<
    java -Xms2000m -jar /usr/picard/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=~{input_bam} \
      OQ=true \
      OUTPUT=~{metrics_filename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    disks: "local-disk " + disk_size + " HDD"
    memory: "3.5 GiB"
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File metrics = "~{metrics_filename}"
  }
}

# Get version of BWA
task GetBwaVersion {
  input {
    Float mem_size_gb = 1
    Int preemptible_tries
    String gotc_docker
    String bwa_path
    String runtime_zones
  }  

  command <<<
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed 
    # because the sed may also fail with that error and that is something we actually want to fail on.
    ~{bwa_path}bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
  >>>

  runtime {
    preemptible: preemptible_tries
    docker: gotc_docker
    memory: "~{mem_size_gb} GiB"
    zones: runtime_zones
  }
  output {
    String version = read_string(stdout())
  }
}

# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
task SamToFastqAndBwaMemAndMba {
  input {
    File input_bam
    String bwa_commandline
    String bwa_version
    String output_bam_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    
    # ref_alt is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
    # listing the reference contigs that are "alternative".
    File ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    Int? mem_size_override
    
    Int compression_level
    Int preemptible_tries
    Int max_retries

    String gotc_docker
    String bwa_path
    String gotc_path
    String runtime_zones
  }
    
  Int mem_size_gb = select_first([mem_size_override, 14])
  Int command_mem_gb = ceil(mem_size_gb/2)

  Float unmapped_bam_size = size(input_bam, "GiB")
  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Float bwa_ref_size = ref_size + size(ref_alt, "GiB") + size(ref_amb, "GiB") + size(ref_ann, "GiB") + size(ref_bwt, "GiB") + size(ref_pac, "GiB") + size(ref_sa, "GiB")
  # Sometimes the output is larger than the input, or a task can spill to disk.
  # In these cases we need to account for the input (1) and the output (1.5) or the input(1), the output(1), and spillage (.5).
  Float disk_multiplier = 2.5
  Int disk_size = ceil(unmapped_bam_size + bwa_ref_size + (disk_multiplier * unmapped_bam_size) + 20)

  command <<<
    set -o pipefail
    set -e

    # set the bash variable needed for the command-line
    bash_ref_fasta=~{ref_fasta}
    # if ref_alt has data in it,
    if [ -s ~{ref_alt} ]; then
      java -Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb}G -jar ~{gotc_path}picard.jar \
        SamToFastq \
        INPUT=~{input_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      ~{bwa_path}~{bwa_commandline} /dev/stdin - 2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) | \
      java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar /usr/gitc/picard.jar \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=~{input_bam} \
        OUTPUT=~{output_bam_basename}.bam \
        REFERENCE_SEQUENCE=~{ref_fasta} \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="~{bwa_version}" \
        PROGRAM_GROUP_COMMAND_LINE="~{bwa_commandline}" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

      grep -m1 "read .* ALT contigs" ~{output_bam_basename}.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    # else ref_alt is empty or could not be found
    else
      exit 1;
    fi
  >>>
  runtime {
    docker: gotc_docker
    preemptible: preemptible_tries
    memory: mem_size_gb + "GB"
    cpu: "16"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: max_retries
    zones: runtime_zones
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File bwa_stderr_log = "~{output_bam_basename}.bwa.stderr.log"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
task SortSam {
  input {
    File input_bam
    String output_bam_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    Int preemptible_tries
    Int compression_level
    Float mem_size_gb = 10

    String gatk_docker
    String gatk_path
    String runtime_zones
  }

  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Float sort_sam_disk_multiplier = 3.25
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20
  
  Int command_mem_gb_sort = ceil(mem_size_gb) - 1
  Int command_mem_gb_fix = ceil((mem_size_gb - 1)/10)

  command <<<
    set -o pipefail

    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb_sort}G" \
      SortSam \
      --INPUT ~{input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
    | \
    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb_fix}G" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ~{output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ~{ref_fasta}
  >>>

  runtime {
    docker: gatk_docker
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "~{mem_size_gb} GB"
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}

# Collect base quality and insert size metrics
task CollectUnsortedReadgroupBamQualityMetrics {
  input {
    File input_bam
    String output_bam_prefix
    Int preemptible_tries
    String runtime_zones
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 20

  command <<<
    java -Xms5000m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=QualityScoreDistribution \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=ALL_READS

    touch ~{output_bam_prefix}.insert_size_metrics
    touch ~{output_bam_prefix}.insert_size_histogram.pdf
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "7 GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File base_distribution_by_cycle_pdf = "~{output_bam_prefix}.base_distribution_by_cycle.pdf"
    File base_distribution_by_cycle_metrics = "~{output_bam_prefix}.base_distribution_by_cycle_metrics"
    File insert_size_histogram_pdf = "~{output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "~{output_bam_prefix}.insert_size_metrics"
    File quality_by_cycle_pdf = "~{output_bam_prefix}.quality_by_cycle.pdf"
    File quality_by_cycle_metrics = "~{output_bam_prefix}.quality_by_cycle_metrics"
    File quality_distribution_pdf = "~{output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "~{output_bam_prefix}.quality_distribution_metrics"
  }
}

# Collect alignment summary and GC bias quality metrics
task CollectReadgroupBamQualityMetrics {
  input {
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int preemptible_tries
    String runtime_zones
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  command <<<
    java -Xms5000m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectGcBiasMetrics \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=READ_GROUP
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "7 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File alignment_summary_metrics = "~{output_bam_prefix}.alignment_summary_metrics"
    File gc_bias_detail_metrics = "~{output_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "~{output_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "~{output_bam_prefix}.gc_bias.summary_metrics"
  }
}

# Collect quality metrics from the aggregated bam
task CollectAggregationMetrics {
  input {
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int preemptible_tries
    String runtime_zones
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  command <<<
    # These are optionally generated, but need to exist for Cromwell's sake
    touch ~{output_bam_prefix}.gc_bias.detail_metrics \
      ~{output_bam_prefix}.gc_bias.pdf \
      ~{output_bam_prefix}.gc_bias.summary_metrics \
      ~{output_bam_prefix}.insert_size_metrics \
      ~{output_bam_prefix}.insert_size_histogram.pdf

    java -Xms5000m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=CollectGcBiasMetrics \
      PROGRAM=QualityScoreDistribution \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    memory: "7 GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
    zones: runtime_zones
  }
  output {
    File alignment_summary_metrics = "~{output_bam_prefix}.alignment_summary_metrics"
    File bait_bias_detail_metrics = "~{output_bam_prefix}.bait_bias_detail_metrics"
    File bait_bias_summary_metrics = "~{output_bam_prefix}.bait_bias_summary_metrics"
    File gc_bias_detail_metrics = "~{output_bam_prefix}.gc_bias.detail_metrics"
    File gc_bias_pdf = "~{output_bam_prefix}.gc_bias.pdf"
    File gc_bias_summary_metrics = "~{output_bam_prefix}.gc_bias.summary_metrics"
    File insert_size_histogram_pdf = "~{output_bam_prefix}.insert_size_histogram.pdf"
    File insert_size_metrics = "~{output_bam_prefix}.insert_size_metrics"
    File pre_adapter_detail_metrics = "~{output_bam_prefix}.pre_adapter_detail_metrics"
    File pre_adapter_summary_metrics = "~{output_bam_prefix}.pre_adapter_summary_metrics"
    File quality_distribution_pdf = "~{output_bam_prefix}.quality_distribution.pdf"
    File quality_distribution_metrics = "~{output_bam_prefix}.quality_distribution_metrics"
  }
}

# Check that the fingerprints of separate readgroups all match
task CrossCheckFingerprints {
  input {
    Array[File] input_bams
    Array[File] input_bam_indexes
    File? haplotype_database_file
    String metrics_filename
    Float total_input_size
    Int preemptible_tries
    String runtime_zones
  }

  Int disk_size = ceil(total_input_size) + 20

  command <<<
    java -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m \
      -jar /usr/picard/picard.jar \
      CrosscheckFingerprints \
      OUTPUT=~{metrics_filename} \
      HAPLOTYPE_MAP=~{haplotype_database_file} \
      EXPECT_ALL_GROUPS_TO_MATCH=true \
      INPUT=~{sep=' INPUT=' input_bams} \
      LOD_THRESHOLD=-20.0
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    preemptible: preemptible_tries
    memory: "2 GB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File metrics = "~{metrics_filename}"
  }
}

# Check that the fingerprint of the sample BAM matches the sample array
task CheckFingerprint {
  input {
    File input_bam
    File input_bam_index
    String output_basename
    File? haplotype_database_file
    File? genotypes
    String sample
    Int preemptible_tries
    String runtime_zones
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 20
  # Picard has different behavior depending on whether or not the OUTPUT parameter ends with a '.', so we are explicitly
  #   passing in where we want the two metrics files to go to avoid any potential confusion.
  String summary_metrics_location = "~{output_basename}.fingerprinting_summary_metrics"
  String detail_metrics_location = "~{output_basename}.fingerprinting_detail_metrics"
  
  command <<<
    java -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms2g  \
      -jar /usr/gitc/picard.jar \
      CheckFingerprint \
      INPUT=~{input_bam} \
      SUMMARY_OUTPUT=~{summary_metrics_location} \
      DETAIL_OUTPUT=~{detail_metrics_location} \
      GENOTYPES=~{genotypes} \
      HAPLOTYPE_MAP=~{haplotype_database_file} \
      SAMPLE_ALIAS="~{sample}" \
      IGNORE_READ_GROUPS=true

  >>>
 runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.1-1540490856"
    preemptible: preemptible_tries
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File summary_metrics = summary_metrics_location
    File detail_metrics = detail_metrics_location
  }
}

# Mark duplicate reads to avoid counting non-independent observations
task MarkDuplicates {
  input {
    Array[File] input_bams
    String output_bam_basename
    String metrics_filename
    Float total_input_size
    Int compression_level
    Int preemptible_tries

    # The program default for READ_NAME_REGEX is appropriate in nearly every case.
    # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
    # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
    String? read_name_regex
    Int memory_multiplier = 1
    Int additional_disk = 20
    
    String gatk_docker
    String gatk_path
    String runtime_zones
  }
  
  # The merged bam will be smaller than the sum of the parts so we need to account for the unmerged inputs and the merged output.
  # Mark Duplicates takes in as input readgroup bams and outputs a slightly smaller aggregated bam. Giving .25 as wiggleroom
  Float md_disk_multiplier = 3
  Int disk_size = ceil(md_disk_multiplier * total_input_size) + additional_disk

  Float memory_size = 7.5 * memory_multiplier
  Int java_memory_size = (ceil(memory_size) - 2)

 # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
 # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
 # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command <<<
    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size}G" \
      MarkDuplicates \
      --INPUT ~{sep=' --INPUT ' input_bams} \
      --OUTPUT ~{output_bam_basename}.bam \
      --METRICS_FILE ~{metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true \
      --CLEAR_DT false \
      --ADD_PG_TAG_TO_READS false
  >>>

  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  input {
    File ref_dict
    Int preemptible_tries
    String python_docker
  }

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("~{ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    docker: python_docker
    preemptible: preemptible_tries
    memory: "2 GB"
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}

task CheckPreValidation {
  input {
    File duplication_metrics
    File chimerism_metrics
    Float max_duplication_in_reasonable_sample
    Float max_chimerism_in_reasonable_sample
    Int preemptible_tries 
    String runtime_zones
  }
   
  command <<<
    set -o pipefail
    set -e

    grep -A 1 PERCENT_DUPLICATION ~{duplication_metrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ~{chimerism_metrics} | grep -v OF_PAIR > chimerism.csv

    python <<CODE

    import csv
    with open('duplication.csv') as dupfile:
      reader = csv.DictReader(dupfile, delimiter='\t')
      for row in reader:
        with open("duplication_value.txt","w") as file:
          file.write(row['PERCENT_DUPLICATION'])
          file.close()

    with open('chimerism.csv') as chimfile:
      reader = csv.DictReader(chimfile, delimiter='\t')
      for row in reader:
        with open("chimerism_value.txt","w") as file:
          file.write(row['PCT_CHIMERAS'])
          file.close()

    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    memory: "2 GB"
    zones: runtime_zones
  }
  output {
    Float duplication_rate = read_float("duplication_value.txt")
    Float chimerism_rate = read_float("chimerism_value.txt")
    Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
  }
}

task ValidateSamFile {
  input {
    File input_bam
    File? input_bam_index
    String report_filename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int? max_output
    Array[String]? ignore
    Boolean? is_outlier_data
    Int? memory_multiplier_override
    Int additional_disk = 20
    Int preemptible_tries
    Int max_retries
    String runtime_zones
  }
  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(1.5 * (size(input_bam, "GiB")) + ref_size) + additional_disk
  
  Int memory_multiplier = select_first([memory_multiplier_override,2])
  Int memory_size = ceil(7 * memory_multiplier)
  Int java_memory_size = (memory_size - 1) * 1000

  command <<<

    java -Xms~{java_memory_size}m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      INPUT=~{input_bam} \
      OUTPUT=~{report_filename} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      ~{"MAX_OUTPUT=" + max_output} \
      IGNORE=~{default="null" sep=" IGNORE=" ignore} \
      MODE=VERBOSE \
      ~{default='SKIP_MATE_VALIDATION=false' true='SKIP_MATE_VALIDATION=true' false='SKIP_MATE_VALIDATION=false' is_outlier_data} \
      IS_BISULFITE_SEQUENCED=false
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    maxRetries: max_retries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File report = "~{report_filename}"
  }
}

# Note these tasks will break if the read lengths in the bam are greater than 250.
task CollectWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    String metrics_filename
    File wgs_coverage_interval_list
    File ref_fasta
    File ref_fasta_index
    Int read_length
    Int preemptible_tries
    String runtime_zones
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  command <<<
    java -Xms2000m -jar /usr/picard/picard.jar \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File metrics = "~{metrics_filename}"
  }
}

# Collect raw WGS metrics (commonly used QC thresholds)
task CollectRawWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    String metrics_filename
    File wgs_coverage_interval_list
    File ref_fasta
    File ref_fasta_index
    Int read_length
    Int preemptible_tries

    Int memory_multiplier = 1
    Int additional_disk = 20  
    String runtime_zones
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + additional_disk

  Int memory_size = ceil((if (disk_size < 110) then 5 else 7) * memory_multiplier)
  String java_memory_size = (memory_size - 1) * 1000
  
  command <<<
    java -Xms~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectRawWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "~{memory_size} GiB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File metrics = "~{metrics_filename}"
  }
}

# Generate a checksum per readgroup
task CalculateReadGroupChecksum {
  input {
    File input_bam
    File input_bam_index
    String read_group_md5_filename
    Int preemptible_tries
    String runtime_zones
  }

  Int disk_size = ceil(size(input_bam, "GiB")) + 20

  command <<<
    java -Xms1000m -jar /usr/picard/picard.jar \
      CalculateReadGroupChecksum \
      INPUT=~{input_bam} \
      OUTPUT=~{read_group_md5_filename}
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "2 GB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File md5_file = "~{read_group_md5_filename}"
  }
}

# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling

task CheckContamination {
  input {
    File input_bam
    File input_bam_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File ref_fasta
    File ref_fasta_index
    String output_prefix
    Int preemptible_tries
    Float contamination_underestimation_factor
    String runtime_zones
  }
  
  Int disk_size = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB")) + 30

  command <<<
    set -e

    # creates a ~{output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ~{output_prefix} \
    --BamFile ~{input_bam} \
    --Reference ~{ref_fasta} \
    --UDPath ~{contamination_sites_ud} \
    --MeanPath ~{contamination_sites_mu} \
    --BedPath ~{contamination_sites_bed} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('~{output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"])/~{contamination_underestimation_factor})
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "7.5 GB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/verify-bam-id:c1cba76e979904eb69c31520a0d7f5be63c72253-1553018888"
    cpu: 2
    zones: runtime_zones
  }
  output {
    File selfSM = "~{output_prefix}.selfSM"
    Float contamination = read_float(stdout())
  }
}

# Convert BAM file to CRAM format
# Note that reading CRAMs directly with Picard is not yet supported
task ConvertToCram {
  input {
    File input_bam
    File ref_fasta
    File ref_fasta_index
    String output_basename
    Int preemptible_tries
    String runtime_zones
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB")
  Int disk_size = ceil(2 * size(input_bam, "GiB") + ref_size) + 20
 
  command <<<
    set -e
    set -o pipefail

    samtools view -C -T ~{ref_fasta} ~{input_bam} | \
    tee ~{output_basename}.cram | \
    md5sum | awk '{print $1}' > ~{output_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    seq_cache_populate.pl -root ./ref/cache ~{ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ~{output_basename}.cram
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    preemptible: preemptible_tries
    memory: "3 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File output_cram = "~{output_basename}.cram"
    File output_cram_index = "~{output_basename}.cram.crai"
    File output_cram_md5 = "~{output_basename}.cram.md5"
  }
}

# Calculates sum of a list of floats
task SumFloats {
  input {
    Array[Float] sizes
    Int preemptible_tries
    String runtime_zones
  }

  command <<<
    python -c "print ~{sep="+" sizes}"
  >>>
  output {
    Float total_size = read_float(stdout())
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    preemptible: preemptible_tries
    zones: runtime_zones
  }
}

task somalier {
  input {
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    File somalier_sites
    String sample_name

    String somalier_docker
    Int preemptible_tries
    String runtime_zones
  }

  Int ref_size = ceil(size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB"))
  Int disk_size = ceil(1.5 * size(input_bam, "GiB") + ref_size + size(somalier_sites,"GiB")) + 20
 
  command <<<
    set -e
    set -o pipefail

    somalier extract --sample-prefix ~{sample_name} --sites ~{somalier_sites} -f ~{ref_fasta} ~{input_bam}
  >>>
  runtime {
    docker: somalier_docker
    preemptible: preemptible_tries
    memory: "3 GB"
    cpu: "1"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File somalier_output  = "~{sample_name}.somalier"
  }
}
