version 1.0
# WORKFLOW DEFINITION 

workflow CallGBAvariants {
  input {
    File cram
    File cram_index
    File ref_fasta
    File ref_fasta_index
    String sample_id
    String? gauch_docker_override
    String gauch_docker = select_first([gauch_docker_override, "zihhuafang/gauchian:v1.0.2"])
    Int? preemptible_attempts_override
    Int preemptible_attempts = select_first([preemptible_attempts_override, "2"])
    String runtime_zones
    Int threads
    }

  call Gauchian {
      input:
        cram = cram,
        cram_index= cram_index,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        sample_id = sample_id,
        gauch_docker = gauch_docker,
	    threads = threads,
        preemptible_attempts = preemptible_attempts,
        runtime_zones = runtime_zones
    }
      # Outputs that will be retained when execution is complete
  output {
      File gauchian_tsv = Gauchian.output_tsv
      File gauchian_json = Gauchian.output_json
  }

}

# TASK DEFINITIONS

task Gauchian {
  input {
  # Command parameters
    File cram
    File cram_index
    File ref_fasta
    File ref_fasta_index
    String sample_id

  # Runtime parameters
    Int additional_disk_space_gb = 10
    Int? machine_mem_gb
    Int preemptible_attempts
    String gauch_docker
    String runtime_zones
    Int threads  
    Int disk_space_gb = ceil(2*size(cram, "G") + size(ref_fasta, "G")) + additional_disk_space_gb
  }

  command {  
  echo "~{cram}" > manifest_gauch.txt
  python -m gauchian --manifest manifest_gauch.txt \
                   --genome 38 \
                   --reference ~{ref_fasta} \
                   --prefix ~{sample_id} \
                   --outDir ~{sample_id} \
                   --threads ~{threads}

  mv ~{sample_id}/*.tsv "~{sample_id}.tsv"
  mv ~{sample_id}/*.json "~{sample_id}.json"
  }
  runtime {
    docker: gauch_docker
    memory: select_first([machine_mem_gb, 10]) + " GB"
    cpu: threads
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    zones: runtime_zones
  }
  output {
    File output_tsv = "~{sample_id}.tsv"
    File output_json = "~{sample_id}.json"
  }
}
