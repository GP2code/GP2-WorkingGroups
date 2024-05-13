version 1.0
#"This pipeline joint-calls GVCFs with GLNexus (https://github.com/dnanexus-rnd/GLnexus). 
# It also permits intervals to be specified so that joint calling only takes place on a subset of intervals"

struct RuntimeAttr {
  Float? mem_gb
  Int? cpu_cores
  Int? disk_gb
  Int? boot_disk_gb
  Int? preemptible_tries
  Int? max_retries
}

workflow GLnexus_jointCall {
  input {
    Array[File] gvcfs
    File dict
    File? bed

    #GLnexus args
    String? config_override
    String config = select_first([config_override, "DeepVariantWGS"])
    Boolean more_PL = false
    Boolean squeeze = false
    Boolean trim_uncalled_alleles = false

    #output args
    String callset_prefix

    #runtime paramenters 
    String? glnexus_docker_override
    String glnexus_docker = select_first([glnexus_docker_override, "quay.io/pacbio/glnexus:v1.4.3"])
    String? runtime_zones_override
    String? gatk_docker_override
    String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:4.2.6.1"])
    String? gatk_path_override
    String gatk_path = select_first([gatk_path_override, "/gatk/gatk"])
    String runtime_zones = select_first([runtime_zones_override,"us-central1-a us-central1-b us-central1-c us-central1-f"])

    RuntimeAttr? runtime_attr_shardgvcf
    RuntimeAttr? runtime_attr_glnexus
    RuntimeAttr? runtime_attr_gather
  }

  parameter_meta {
    gvcfs:       "gVCF files to perform joint calling upon"
    dict:        "reference sequence dictionary"
    bed:         "intervals to which joint calling should be restricted"
    callset_prefix: "output prefix for joined-called VCF files"

    config:      "configuration preset name or .yml filename"
    more_PL:     "include PL from reference bands and other cases omitted by default"
    squeeze:     "reduce pVCF size by suppressing detail in cells derived from reference bands"
    trim_uncalled_alleles: "remove alleles with no output GT calls in postprocessing"

    runtime_attr_override: "override default runtime attributes"
  }
    
  # List all of the contigs in the reference
  call GetRanges { 
    input: 
      dict = dict, 
      bed = bed,
      glnexus_docker = glnexus_docker,
      runtime_zones = runtime_zones
  }

  # Shard all gVCFs into intervals
  scatter (p in gvcfs) {
    call ShardVCFByRanges { 
      input: 
        gvcf = p, 
        tbi = p + ".tbi", 
        ranges = GetRanges.ranges,
        docker = glnexus_docker,
        runtime_zones = runtime_zones,
        runtime_attr_override =runtime_attr_shardgvcf
    }
  }

  # Joint-call in parallel over chromosomes
  scatter (idx in range(length(ShardVCFByRanges.sharded_gvcfs[0]))) {
    Array[File] per_contig_gvcfs = transpose(ShardVCFByRanges.sharded_gvcfs)[idx]

    call glnexus_call {
      input:
        gvcfs = per_contig_gvcfs,
        config = config,
        more_PL = more_PL,
        squeeze = squeeze,
        trim_uncalled_alleles = trim_uncalled_alleles,
        prefix = callset_prefix + ".raw.",
        docker = glnexus_docker,
        runtime_zones = runtime_zones,
        runtime_attr_override =runtime_attr_glnexus
      }
  }

  output {
    Array[File] joint_vcf = glnexus_call.joint_vcf
  }
}

#Task definition
#Generate interval for shard vcf
task GetRanges {
  input {
      File dict
      File? bed
      Int preemptible_tries  
      String glnexus_docker
      String runtime_zones
   }
 
  Int disk_size = 1 + ceil(size(dict, "GB"))

  command <<<
        set -euxo pipefail

        if [[ "~{defined(bed)}" == "true" ]]; then
            cat ~{bed} | awk '{ print $1 ":" $2 "-" $3 }' > ranges.txt
        else
            grep '^@SQ' ~{dict} | \
                awk '{ print $2, $3 }' | \
                sed 's/[SL]N://g' | \
                awk '{ print $1 ":0-" $2 }' | \
                awk 'FNR < 25' \
                > ranges.txt
        fi
  >>>

  output {
    Array[String] ranges = read_lines("ranges.txt")
  }
    
  runtime {
    memory: "3 GiB"
    docker: glnexus_docker
    zones: runtime_zones
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

#Joint-call gVCFs with GLNexus
task glnexus_call {
  input {
    Array[File] gvcfs

    String config
    Boolean more_PL = false
    Boolean squeeze = false
    Boolean trim_uncalled_alleles = false

    Int? num_cpus_override
    # for a small number of shards 7.5 GB (n1-standard-2) was not enough, so we're increasing to 13 GB (n1-highmem-2)
    Float mem_default = 7.5
    String prefix
    String runtime_zones
    String docker
    Int? additional_disk_override 

    RuntimeAttr? runtime_attr_override
  }
    Int additional_disk = select_first([additional_disk_override, 20])
    Int num_cpus = select_first([num_cpus_override, 2])
    Int disk_size = additional_disk + 5*ceil(size(gvcfs, "GiB"))
    Float mem = if defined(num_cpus_override) then 3.75*num_cpus else mem_default
    
  command <<<
      set -x

      # For guidance on performance settings, see https://github.com/dnanexus-rnd/GLnexus/wiki/Performance
      ulimit -Sn 65536

      echo ~{gvcfs[0]} | sed 's/.*locus_//' | sed 's/.g.vcf.gz//' | sed 's/___/\t/g' > range.bed
      chr=$(echo ~{gvcfs[0]} | sed 's/.*locus_//' | sed 's/.g.vcf.gz//' | sed 's/___/\t/g' | awk '{print $1}') 

      glnexus_cli \
          --config ~{config} \
          --bed range.bed \
          ~{if more_PL then "--more-PL" else ""} \
          ~{if squeeze then "--squeeze" else ""} \
          ~{if trim_uncalled_alleles then "--trim-uncalled-alleles" else ""} \
          --list ~{write_lines(gvcfs)} \
       > ~{prefix}${chr}.bcf
      
  >>>

  output {
      File joint_vcf = glob("~{prefix}*.bcf")
  }

  #########################
  RuntimeAttr default_attr = object {
      cpu_cores:          num_cpus,
      mem_gb:             mem,
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
      cpu: runtime_attr.cpu_cores
      memory: runtime_attr.mem_gb + " GiB"
      disks: "local-disk " +  runtime_attr.disk_gb + " SSD"
      bootDiskSizeGb: runtime_attr.boot_disk_gb
      preemptible: runtime_attr.preemptible_tries
      maxRetries: runtime_attr.max_retries
      zones: runtime_zones
  }
}