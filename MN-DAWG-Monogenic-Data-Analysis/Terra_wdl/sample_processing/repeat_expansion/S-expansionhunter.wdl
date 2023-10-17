version 1.0
# WDL for calling target repeats by Expansion Hunter (https://github.com/Illumina/ExpansionHunter)
# Inputs: bam/cram
#         sample_id
#         sex
#         varaint catalog: https://github.com/Illumina/ExpansionHunter/tree/master/variant_catalog/hg38

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

workflow ExpansionHunter {
    input {
        File bams_or_cram
        File bams_or_cram_index
        String sample_id
        String sex
        File reference_fasta
        File reference_fasta_index
        File variant_catalog_json
        String expansion_hunter_docker
        String python_docker

        RuntimeAttr? runtime_eh
        RuntimeAttr? runtime_concat
    }

    parameter_meta {
        sexes: "sex of the sample (in dataable)"
        sample_ids: "One ID per sample, in the same order as the files in bams_or_crams. These IDs must match the ID given in the second column (`Individual ID` column) of the given PED file. These IDs will also be used as an output prefix."
    }

    call RunExpansionHunter {
      input:
        bam_or_cram = bams_or_cram,
        bam_or_cram_index = bams_or_cram_index,
        reference_fasta = reference_fasta,
        reference_fasta_index = reference_fasta_index,
        sample_id = sample_id,
        sex = sex,
        variant_catalog = variant_catalog_json,
        expansion_hunter_docker = expansion_hunter_docker,
        runtime_override = runtime_eh
    }


    output {
        File variants_tsv = RunExpansionHunter.variants_tsv
        File alleles_tsv = RunExpansionHunter.alleles_tsv
        File vcf_gz = RunExpansionHunter.vcf_gz
        File vcf_gz_index = RunExpansionHunter.vcf_gz_index
        File overlapping_reads = RunExpansionHunter.overlapping_reads
        File overlapping_reads_index = RunExpansionHunter.overlapping_reads_index
    }
}


task RunExpansionHunter {
    input {
        File bam_or_cram
        File bam_or_cram_index
        File reference_fasta
        File reference_fasta_index
        File variant_catalog
        String sample_id
        String sex
        File? ped_file
        String expansion_hunter_docker
        RuntimeAttr? runtime_override
    }

    command <<<
        set -euxo pipefail

        REF="$(basename "~{reference_fasta}")"
        mv ~{reference_fasta} $REF
        mv ~{reference_fasta_index} $REF.fai

        
        ExpansionHunter \
            --reads ~{bam_or_cram} \
            --reference $REF \
            --variant-catalog ~{variant_catalog} \
            --output-prefix ~{sample_id} \
            --cache-mates \
            --sex ~{sex}

        #sort bam
        samtools sort -o ~{sample_id}_realigned_sorted.bam ~{sample_id}_realigned.bam
        samtools index ~{sample_id}_realigned_sorted.bam

        bgzip ~{sample_id}.vcf && tabix -p vcf ~{sample_id}.vcf.gz

        python /opt/str/combine_expansion_hunter_json_to_tsv.py -o ~{sample_id} ~{sample_id}.json
        mv ~{sample_id}.*_json_files_alleles.tsv ~{sample_id}_alleles.tsv
        mv ~{sample_id}.*_json_files_variants.tsv ~{sample_id}_variants.tsv
    >>>

    RuntimeAttr runtime_default = object {
        cpu_cores: 1,
        mem_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 1,
        disk_gb: 10 + ceil(size([
            bam_or_cram,
            bam_or_cram_index,
            reference_fasta,
            reference_fasta_index], "GiB"))
    }
    RuntimeAttr runtime_attr = select_first([runtime_override, runtime_default])

    runtime {
        docker: expansion_hunter_docker
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb])  + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    }

    output {
        File variants_tsv = "${sample_id}_variants.tsv"
        File alleles_tsv = "${sample_id}_alleles.tsv"
        File vcf_gz = "${sample_id}.vcf.gz"
        File vcf_gz_index = "${sample_id}.vcf.gz.tbi"
        File overlapping_reads = "${sample_id}_realigned_sorted.bam"
        File overlapping_reads_index = "${sample_id}_realigned_sorted.bam.bai"
    }
}