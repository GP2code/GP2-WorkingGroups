## Copyright Broad Institute, 2018
## 
## This WDL converts paired FASTQ to uBAM and adds read group information 
##
## Requirements/expectations :
## - Pair-end sequencing data in FASTQ format (one file per orientation)
## - One or more read groups, one per pair of FASTQ files  
## - A readgroup.list file with the following format :  
##   ``readgroup   fastq_pair1   fastq_pair2   sample_name   library_name   platform_unit   run_date   platform_name   sequecing_center``
##
## Outputs :
## - Set of unmapped BAMs, one per read group
## - File of a list of the generated unmapped BAMs
##
## Cromwell version support 
## - Successfully tested on v32
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# WORKFLOW DEFINITION
workflow ConvertPairedFastQsToUnmappedBamWf {
  File readgroup_list  
  Array[Array[String]] readgroup_array = read_tsv(readgroup_list) 
  String ubam_list_name = basename(readgroup_list,".txt") + "_ubams"
  
  String? gatk_docker_override
  String gatk_docker = select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.2.6.0"])
  String? gatk_path_override
  String gatk_path = select_first([gatk_path_override, "gatk"])
  Int? preemptible_attempts
  String runtime_zones

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
        preemptible_attempts = preemptible_attempts,
        runtime_zones = runtime_zones
    }
  }

  #Create a file with a list of the generated ubams
  call CreateFoFN {
    input:
      array_of_files = PairedFastQsToUnmappedBAM.output_bam,
      fofn_name = ubam_list_name
  }

  # Outputs that will be retained when execution is complete
  output {
    Array[File] output_bams = PairedFastQsToUnmappedBAM.output_bam
    File unmapped_bam_list = CreateFoFN.fofn_list
  }
}

# TASK DEFINITIONS

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
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
  Int additional_disk_space_gb= select_first([additional_disk_gb_override, 50])
  Int? machine_mem_gb
  Int? preemptible_attempts
  String gatk_docker
  String gatk_path
  String runtime_zones
  
  Int disk_space_gb = ceil(size(fastq_1, "GB")* 2 + size(fastq_2, "GB") * 2) + additional_disk_space_gb
  command {  
    ${gatk_path} --java-options "-Xmx3000m" \
    FastqToSam \
    --FASTQ ${fastq_1} \
    --FASTQ2 ${fastq_2} \
    --OUTPUT ${readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ${readgroup_name} \
    --SAMPLE_NAME ${sample_name} \
    --LIBRARY_NAME ${library_name} \
    --PLATFORM_UNIT ${platform_unit} \
    --RUN_DATE ${run_date} \
    --PLATFORM ${platform_name} \
    --SEQUENCING_CENTER ${sequencing_center} 
  }
  runtime {
    docker: gatk_docker
    memory: select_first([machine_mem_gb, 10]) + " GB"
    cpu: "1"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: select_first([preemptible_attempts, 3])
    zones: runtime_zones
  }
  output {
    File output_bam = "${readgroup_name}.unmapped.bam"
  }
}

task CreateFoFN {
  # Command parameters
  Array[String] array_of_files
  String fofn_name
  
  command {
    mv ${write_lines(array_of_files)}  ${fofn_name}.list
  }
  output {
    File fofn_list = "${fofn_name}.list"
  }
  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }
}