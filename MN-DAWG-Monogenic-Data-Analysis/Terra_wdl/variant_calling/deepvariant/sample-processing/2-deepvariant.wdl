version 1.0

struct RuntimeAttr {
  Float? mem_gb
  Int? cpu_cores
  Int? disk_gb
  Int? boot_disk_gb
  Int? preemptible_tries
  Int? max_retries
}

workflow DeepVariant {
    input {
        File referenceFasta
        File referenceFastaIndex
        File ref_dict
        File input_cram_or_bam
        File? input_cram_or_bam_index
        String sample_name
        String sex

        #deepvariant args
        String model_type
        File PAR_region

        # gvcf/vcf QC files 
        File dbSNP_vcf
        File dbSNP_vcf_index
        Array[File] known_indels_sites_VCFs
        Array[File] known_indels_sites_indices
        File wgs_calling_interval_list
        File wgs_evaluation_interval_list
        
        # inputs for sharding gvcf
        File dict
        File? bed
        File unpadded_intervals_file 
        File regions = select_first([bed,unpadded_intervals_file])
        
        #runtime paramters
        RuntimeAttr? runtime_attr_deepvariant
        RuntimeAttr? runtime_attr_shardgvcf
        String? r_zones_override
        String r_zones = select_first([r_zones_override, "us-central1-a us-central1-b us-central1-c us-central1-f"])
        String? deepvariant_docker_override
        # deepvariant 1.6.0 
        String deepvariant_docker = select_first([deepvariant_docker_override, "google/deepvariant:CL590726281"])
        String? gatk_docker_override
        String gatk_docker = select_first([gatk_docker_override, "us.gcr.io/broad-gatk/gatk:4.2.6.0"])
        String? gatk_path_override
        String gatk_path = select_first([gatk_path_override, "gatk"])
        String? glnexus_docker_override
        String glnexus_docker = select_first([glnexus_docker_override, "quay.io/pacbio/glnexus:v1.4.3"])
        Int preemptible_tries = 3
    }
        
    Boolean is_bam =
      basename(input_cram_or_bam, ".bam") + ".bam" ==
      basename(input_cram_or_bam)

    File bam_or_cram_index =
      if defined(input_cram_or_bam_index) then
        select_first([input_cram_or_bam_index])
      else
        input_cram_or_bam + if is_bam then ".bai" else ".crai"
        
    call deepvariant {
      input:
        input_cram_or_bam = input_cram_or_bam,
        input_cram_or_bam_index = bam_or_cram_index,
        sample_name = sample_name,
        sex = sex,
        model_type = model_type,
        PAR_region = PAR_region,
        referenceFasta = referenceFasta,
        referenceFastaIndex = referenceFastaIndex,
        deepvariant_docker = deepvariant_docker,
        runtime_attr_override = runtime_attr_deepvariant,
        r_zones = r_zones
    }

    # QC the GVCF
    call CollectGvcfCallingMetrics {
      input:
        input_vcf = deepvariant.outputGVCF,
        input_vcf_index = deepvariant.outputGVCFIndex,
        metrics_basename = sample_name,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        ref_dict = ref_dict,
        wgs_evaluation_interval_list = wgs_evaluation_interval_list,
        preemptible_tries = preemptible_tries,
        runtime_zones = r_zones
    }

    Int num_of_original_intervals = length(read_lines(regions))    
    Int merge_count = 3
    
    call DynamicallyCombineIntervals {
      input:
        intervals = unpadded_intervals_file,
        merge_count = merge_count,
        runtime_zones = r_zones
    }

    Array[String] ranges = read_lines(DynamicallyCombineIntervals.output_intervals)

    # Shard all gVCFs into intervals
    call ShardVCFByRanges { 
      input: 
        gvcf = deepvariant.outputGVCF, 
        tbi = deepvariant.outputGVCFIndex, 
        ranges = ranges,
        docker = glnexus_docker,
        runtime_zones = r_zones,
        runtime_attr_override =runtime_attr_shardgvcf
    }

    call Get_gvcf_fofn {
      input:
        input_gvcfs = ShardVCFByRanges.sharded_gvcfs,
        fofn_name = sample_name +"_sharded_gvcfs.list"
    }
   
    output {
        File sample_gvcf=deepvariant.outputGVCF
        File sample_vcf=deepvariant.outputVCF
        File gvcf_summary_metrics = CollectGvcfCallingMetrics.summary_metrics
        File gvcf_detail_metrics = CollectGvcfCallingMetrics.detail_metrics
        File gvcf_fofn = Get_gvcf_fofn.fofn_list
    }
}

task deepvariant {
    input {
        File referenceFasta
        File referenceFastaIndex
        File input_cram_or_bam
        File input_cram_or_bam_index
        String model_type
        String sample_name
        String sex
        File PAR_region
        Boolean? VCFStatsReport_override

        # deepvariant optional arguments
        Array[String]? make_exampleExtraArgs
        Array[String]? callvariantsExtraArgs
        Array[String]? postprocessVariantsExtraArgs
        
        #runtime paramters
        String r_zones
        String deepvariant_docker
        RuntimeAttr? runtime_attr_override
        Int? num_cpus_override
        Int? mem_override
        Int? additional_disk_override
    }
    # output
    String outputGVcf = sample_name + ".g.vcf.gz"
    String outputVcf =  sample_name + ".vcf.gz"
    
    Boolean is_male = sex == "male"
    Boolean VCFStatsReport = select_first([VCFStatsReport_override, false])
    # other variables (n1-standard-16)
    Int cram_size = round(size(input_cram_or_bam, "GiB")) 

    Int num_cpus = select_first([num_cpus_override, 16])
    Int machine_mem = select_first([mem_override, 80])

    Int additional_disk = select_first([additional_disk_override, 40])
    Int disk_size = round(2 * size(input_cram_or_bam, "GiB")) + additional_disk
 
    output {
        File outputVCF = outputVcf
        File outputVCFIndex = outputVcf + ".tbi"
        File outputGVCF = outputGVcf
        File outputGVCFIndex = outputGVcf + ".tbi"
        Array[File]? outputVCFStatsReport = glob("*.visual_report.html")
    }

    command <<<
      set -e
      make_example=$(echo ~{sep=',' make_exampleExtraArgs})
      callvariants=$(echo ~{sep=',' callvariantsExtraArgs})
      postprocessVariants=$(echo ~{sep=',' postprocessVariantsExtraArgs})
      haploid=$(echo "\"chrX,chrY\"")
      PAR_bed=$(echo "~{PAR_region}")
      
      /opt/deepvariant/bin/run_deepvariant \
        --ref=~{referenceFasta} \
        --reads=~{input_cram_or_bam} \
        --model_type=~{model_type} \
        --output_vcf=~{outputVcf} \
        --output_gvcf=~{outputGVcf} \
        --intermediate_results_dir=tmp \
        --num_shards=$(nproc) \
        --sample_name=~{sample_name} \
        ~{if defined(make_exampleExtraArgs) then "--make_examples_extra_args $make_example" else ""} \
        ~{if defined(callvariantsExtraArgs) then "--call_variants_extra_args $callvariants" else ""} \
        ~{if defined(postprocessVariantsExtraArgs) then "--postprocess_variants_extra_args $postprocessVariants" else ""} \
        ~{true="--vcf_stats_report" false="--novcf_stats_report" VCFStatsReport} \
        ~{true="--haploid_contigs=$haploid" false="" is_male} \
        ~{true="--par_regions_bed=$PAR_bed" false="" is_male}
    >>>

    RuntimeAttr runtime_attr_default = object {
      cpu_cores: num_cpus,
      mem_gb: machine_mem,
      boot_disk_gb: 20,
      preemptible_tries: 5,
      max_retries: 1,
      disk_gb: disk_size
    }
    
    RuntimeAttr runtime_attr = select_first([
      runtime_attr_override,
      runtime_attr_default])

    runtime {
      docker: deepvariant_docker
      cpu: runtime_attr.cpu_cores
      memory: runtime_attr.mem_gb + " GiB"
      disks: "local-disk " + runtime_attr.disk_gb + " SSD"
      bootDiskSizeGb: runtime_attr.boot_disk_gb
      preemptible: runtime_attr.preemptible_tries
      maxRetries: runtime_attr.max_retries
      zones: r_zones
    }
}

# Collect variant calling metrics from GVCF output
task CollectGvcfCallingMetrics {
  input {
    File input_vcf
    File input_vcf_index
    String metrics_basename
    File dbSNP_vcf
    File dbSNP_vcf_index
    File ref_dict
    File wgs_evaluation_interval_list
    Int preemptible_tries
    String runtime_zones
  }

  Int disk_size = ceil(size(input_vcf, "GiB") + size(dbSNP_vcf, "GiB")) + 20

  command <<<
    java -Xms2000m -jar /usr/picard/picard.jar \
      CollectVariantCallingMetrics \
      INPUT=~{input_vcf} \
      OUTPUT=~{metrics_basename} \
      DBSNP=~{dbSNP_vcf} \
      SEQUENCE_DICTIONARY=~{ref_dict} \
      TARGET_INTERVALS=~{wgs_evaluation_interval_list} \
      GVCF_INPUT=true
  >>>

  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8"
    preemptible: preemptible_tries
    memory: "3 GB"
    disks: "local-disk " + disk_size + " HDD"
    zones: runtime_zones
  }
  output {
    File summary_metrics = "~{metrics_basename}.variant_calling_summary_metrics"
    File detail_metrics = "~{metrics_basename}.variant_calling_detail_metrics"
  }
}

task DynamicallyCombineIntervals {
  input {
    File intervals
    Int merge_count
    Int max_retries = 1
    Int preemptible_tries = 3 
    String runtime_zones
  }
  
  command <<<
    python << CODE
    def parse_interval(interval):
        colon_split = interval.split(":")
        chromosome = colon_split[0]
        dash_split = colon_split[1].split("-")
        start = int(dash_split[0])
        end = int(dash_split[1])
        return chromosome, start, end

    def add_interval(chr, start, end):
        lines_to_write.append(chr + ":" + str(start) + "-" + str(end))
        return chr, start, end

    count = 0
    chain_count = ~{merge_count}
    l_chr, l_start, l_end = "", 0, 0
    lines_to_write = []
    with open("~{intervals}") as f:
        with open("out.intervals", "w") as f1:
            for line in f.readlines():
                # initialization
                if count == 0:
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue
                # reached number to combine, so spit out and start over
                if count == chain_count:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
                    continue

                c_chr, c_start, c_end = parse_interval(line)
                # if adjacent keep the chain going
                if c_chr == w_chr and c_start == w_end + 1:
                    w_end = c_end
                    count += 1
                    continue
                # not adjacent, end here and start a new chain
                else:
                    l_char, l_start, l_end = add_interval(w_chr, w_start, w_end)
                    w_chr, w_start, w_end = parse_interval(line)
                    count = 1
            if l_char != w_chr or l_start != w_start or l_end != w_end:
                add_interval(w_chr, w_start, w_end)
            f1.writelines("\n".join(lines_to_write))
    CODE
    >>>

  runtime {
    memory: "3 GiB"
    maxRetries: max_retries
    preemptible: preemptible_tries
    docker: "python:2.7"
    zones: runtime_zones
  }

  output {
    File output_intervals = "out.intervals"
  }
}

#Split VCF into smaller ranges for parallelization
task ShardVCFByRanges {
  input {
    File gvcf
    File tbi
    Array[String] ranges

    String docker
    String runtime_zones

    RuntimeAttr? runtime_attr_override
  }
    String gvcf_basename= basename(gvcf, ".g.vcf.gz")
    Int disk_size = 10 + 2*ceil(size(gvcf, "GiB"))

  command <<<
      set -oe pipefail

      mkdir per_interval

      INDEX=0
      for RANGE in ~{sep=' ' ranges}
      do
          PINDEX=$(printf "%06d" $INDEX)
          FRANGE=$(echo $RANGE | sed 's/[:-]/___/g')
          OUTFILE="per_interval/$PINDEX.~{gvcf_basename}.locus_$FRANGE.g.vcf.gz"


          bcftools view ~{gvcf} $RANGE | bgzip > $OUTFILE

          INDEX=$(($INDEX+1))
      done
  >>>

  output {
      Array[File] sharded_gvcfs = glob("per_interval/*.g.vcf.gz")
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          1,
      mem_gb:             3.75,
      disk_gb:            disk_size,
      boot_disk_gb:       10,
      preemptible_tries:  3,
      max_retries:        1
  }

  RuntimeAttr runtime_attr = select_first([
      runtime_attr_override,
      default_attr])

  runtime {
      docker: docker
      cpu: 1
      memory: 1 + " GiB"
      disks: "local-disk " +  runtime_attr.disk_gb + " SSD"
      bootDiskSizeGb: runtime_attr.boot_disk_gb
      preemptible: runtime_attr.preemptible_tries
      maxRetries: runtime_attr.max_retries
      zones: runtime_zones
  }
}

task Get_gvcf_fofn {
  input {
    # Command parameters
    Array[String] input_gvcfs
    String fofn_name
  }
  
  command <<<
     mv ~{write_lines(input_gvcfs)}  ~{fofn_name}
  >>>
  output {
    File fofn_list = "~{fofn_name}"
  }
  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }
}