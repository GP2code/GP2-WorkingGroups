version 1.0
## Copyright Broad Institute, 2019
## Purpose: 
## Generate a sample_map file, which can be used for JointGenotyping workflow
##
## Requirements/expectations :
## - An array of file paths
## - An array of file names
## - Name of output sample_map
##
## Outputs :
## - sample map file
##
## Cromwell version support 
## - Successfully tested on v47
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.



# WORKFLOW DEFINITION 

workflow GenerateSampleMap {
  input {
    Array[String] sample_names
    Array[String] file_paths
    String sample_map_name
  }
  
  call GenerateSampleMapFile {
    input:
      sample_names = sample_names,
      file_paths = file_paths,
      outfile = sample_map_name + ".sample_map"
  }
  
  output {
    File sample_map = GenerateSampleMapFile.sample_map
  }

}


# TASK DEFINITIONS
task GenerateSampleMapFile {
  input{
    # Command parameters
    Array[String] sample_names
    Array[String] file_paths
    String outfile 
    
    # Runtime parameters
    String docker = "python:latest"
    Int machine_mem_gb = 7
    Int disk_space_gb = 100
    Int preemptible_attempts = 3
  }
    command <<<
    set -oe pipefail
    
    python << CODE
    file_paths = ['~{sep="','" file_paths}']
    sample_names = ['~{sep="','" sample_names}']

    if len(file_paths)!= len(sample_names):
      print("Number of File Paths does not Equal Number of File Names")
      exit(1)

    with open("sample_map_file.txt", "w") as fi:
      for i in range(len(file_paths)):
        fi.write(sample_names[i] + "\t" + file_paths[i] + "\n") 

    CODE
    mv sample_map_file.txt ~{outfile}
    >>>

    runtime {
        docker: docker
        memory: machine_mem_gb + " GB"
        disks: "local-disk " + disk_space_gb + " HDD"
        preemptible: 3
    }

    output {
        File sample_map = outfile
    }
}