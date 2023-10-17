version 1.0

###WDL version of parliment2 for single sample SV calling 
#  https://github.com/slzarate/parliament2

workflow Parliament2 {
    input {
        File alignmentFile
        File alignmentIndex
        File refFasta 
        File refIndex
        Boolean filterContigs
        Boolean runBreakdancer
        Boolean runBreakseq
        Boolean runCNVnator
        Boolean runDelly
        Boolean runLumpy
        Boolean runManta
        Boolean genotypeVCFs
        Boolean runSURVIVOR
        Boolean runJasmine
        Boolean visualizeSVs
        ## Runtime parameters
        String runtime_zones 
        Int max_retries = 2
        Int preemptible_tries = 3
    }

    call P2Prep {
        input:
            inputFile = alignmentFile,
            inputIndex = alignmentIndex,
            refFasta = refFasta,
            filterContigs = filterContigs,
            runtime_zones = runtime_zones,
            preemptible_tries = preemptible_tries
    }

    Array[String] chromosomes = read_lines(P2Prep.contigs)

    if (runBreakdancer) {
        call PrepareBreakdancer {
            input:
                inputBam = P2Prep.outputBam,
                runtime_zones = runtime_zones,
                preemptible_tries = preemptible_tries
        }

        scatter (chromosome in chromosomes) {
            call BreakdancerChrom {
                input:
                    inputBam = P2Prep.outputBam,
                    inputBai = P2Prep.outputBai,
                    configFile = PrepareBreakdancer.configFile,
                    contig = chromosome,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }

        call GatherBreakdancer {
            input:
                breakdancerCtx = BreakdancerChrom.breakdancerOutput,
                bamBase = P2Prep.bamBase,
                runtime_zones = runtime_zones,
                preemptible_tries = preemptible_tries
        }

        if (genotypeVCFs) {
            call SVTyper as typeBreakdancer {
                input:
                    inputBam = P2Prep.outputBam,
                    inputBai = P2Prep.outputBai,
                    refFasta = refFasta,
                    refIndex = refIndex,
                    inputVCF = GatherBreakdancer.breakdancerVCF,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }
    }
        
    if (runCNVnator) {
        scatter (chromosome in chromosomes) {
            call CNVnatorChrom {
                input:
                    inputBam = P2Prep.outputBam,
                    inputBai = P2Prep.outputBai,
                    refFasta = refFasta,
                    refIndex = refIndex,
                    contig = chromosome,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }

        call GatherCNVnator {
            input:
                CNVnatorCalls = CNVnatorChrom.CNVnatorOutput,
                bamBase = P2Prep.bamBase,
                runtime_zones = runtime_zones,
                preemptible_tries = preemptible_tries
            }

        if (genotypeVCFs) {
            call SVTyper as typeCNVnator {
                input:
                    inputBam = P2Prep.outputBam,
                    inputBai = P2Prep.outputBai,
                    refFasta = refFasta,
                    refIndex = refIndex,
                    inputVCF = GatherCNVnator.CNVnatorVCF,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }
    }

    if (runDelly || runLumpy) {
        scatter (chromosome in chromosomes) {
            call SplitByChrom {
                input:
                    inputBam = P2Prep.outputBam,
                    inputBai = P2Prep.outputBai,
                    contig = chromosome,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }

        Array[Pair[File, File]] indexedBams = zip(SplitByChrom.splitBam, SplitByChrom.splitBai)

        if (runDelly) {
            scatter (delly_i in indexedBams) {
                call DellyChrom {
                    input:
                        inputBam = delly_i.left,
                        inputBai = delly_i.right,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        runtime_zones = runtime_zones,
                        max_retries = max_retries,
                        preemptible_tries = preemptible_tries
                }
            }
            call GatherDelly {
                input:
                    dellyDelVCFs = DellyChrom.dellyDeletion,
                    dellyDupVCFs = DellyChrom.dellyDuplication,
                    dellyInsVCFs = DellyChrom.dellyInsertion,
                    dellyInvVCFs = DellyChrom.dellyInversion,
                    bamBase = P2Prep.bamBase,
                    runtime_zones = runtime_zones,
                    preemptible_tries = preemptible_tries
            }

            if (genotypeVCFs) {
                call SVTyper as typeDellyDel {
                    input:
                        inputBam = P2Prep.outputBam,
                        inputBai = P2Prep.outputBai,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        inputVCF = GatherDelly.dellyDeletion,
                        runtime_zones = runtime_zones,
                        max_retries = max_retries,
                        preemptible_tries = preemptible_tries
                    }

                call SVTyper as typeDellyDup {
                    input:
                        inputBam = P2Prep.outputBam,
                        inputBai = P2Prep.outputBai,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        inputVCF = GatherDelly.dellyDuplication,
                        runtime_zones = runtime_zones,
                        max_retries = max_retries,
                        preemptible_tries = preemptible_tries
                }

                call SVTyper as typeDellyIns {
                    input:
                        inputBam = P2Prep.outputBam,
                        inputBai = P2Prep.outputBai,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        inputVCF = GatherDelly.dellyInsertion,
                        runtime_zones = runtime_zones,
                        max_retries = max_retries,
                        preemptible_tries = preemptible_tries
                }

                call SVTyper as typeDellyInv {
                    input:
                        inputBam = P2Prep.outputBam,
                        inputBai = P2Prep.outputBai,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        inputVCF = GatherDelly.dellyInversion,
                        runtime_zones = runtime_zones,
                        max_retries = max_retries,
                        preemptible_tries = preemptible_tries
                }
            }
        }

        if (runLumpy) {
            scatter (lumpy_i in indexedBams) {
                call LumpyChrom {
                    input:
                        inputBam = lumpy_i.left,
                        inputBai = lumpy_i.right,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        runtime_zones = runtime_zones,
                        preemptible_tries = preemptible_tries
                }
            }

            call GatherLumpy {
                input:
                    lumpyVCFs = LumpyChrom.lumpyOutput,
                    bamBase = P2Prep.bamBase,
                    runtime_zones = runtime_zones,
                    preemptible_tries = preemptible_tries
            }

            if (genotypeVCFs) {
                call SVTyper as typeLumpy {
                    input:
                        inputBam = P2Prep.outputBam,
                        inputBai = P2Prep.outputBai,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        inputVCF = GatherLumpy.lumpyVCF,
                        runtime_zones = runtime_zones,
                        max_retries = max_retries,
                        preemptible_tries = preemptible_tries
                }
            }
        }
    }

    if (runManta) {
        call Manta {
            input:
                inputBam = P2Prep.outputBam,
                inputBai = P2Prep.outputBai,
                refFasta = refFasta,
                refIndex = refIndex,
                contigs = P2Prep.contigs,
                runtime_zones = runtime_zones,
                max_retries = max_retries,
                preemptible_tries = preemptible_tries
        }

        if (genotypeVCFs) {
            call SVTyper as typeManta {
                input:
                    inputBam = P2Prep.outputBam,
                    inputBai = P2Prep.outputBai,
                    refFasta = refFasta,
                    refIndex = refIndex,
                    inputVCF = Manta.mantaVCF,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }
    }

    if (runBreakseq) {
        call Breakseq { 
            input: 
                inputBam = P2Prep.outputBam,
                inputBai = P2Prep.outputBai,
                refFasta = refFasta,
                refIndex = refIndex,
                reference = P2Prep.refType,
                runtime_zones = runtime_zones,
                max_retries = max_retries,
                preemptible_tries = preemptible_tries
        }

        if (genotypeVCFs) {
            call SVTyper as typeBreakseq {
                input:
                    inputBam = P2Prep.outputBam,
                    inputBai = P2Prep.outputBai,
                    refFasta = refFasta,
                    refIndex = refIndex,
                    inputVCF = Breakseq.breakseqVCF,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
                }
        }
    }

    if (genotypeVCFs) {
        Array[File?] potentialSVtyped = [ typeBreakdancer.svtypedVCF,
                                          typeBreakseq.svtypedVCF,
                                          typeCNVnator.svtypedVCF,
                                          typeDellyDel.svtypedVCF,
                                          typeDellyDup.svtypedVCF,
                                          typeDellyIns.svtypedVCF,
                                          typeDellyInv.svtypedVCF,
                                          typeLumpy.svtypedVCF,
                                          typeManta.svtypedVCF
                                        ]
        
        Array[File] svtypedVCFs = select_all(potentialSVtyped)

        if (runSURVIVOR) {
            call SURVIVOR as svtypedSurvivor {
                input:
                    vcfs = svtypedVCFs,
                    bamBase = P2Prep.bamBase,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }

        if (runJasmine) {
            call Jasmine as svtypedJasmine {
                input:
                    vcfs = svtypedVCFs,
                    bamBase = P2Prep.bamBase,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }
    }

    if (!genotypeVCFs) {
        Array[File?] potentialVCF = [ GatherBreakdancer.breakdancerVCF,
                                      Breakseq.breakseqVCF,
                                      GatherCNVnator.CNVnatorVCF,
                                      GatherDelly.dellyDeletion,
                                      GatherDelly.dellyDuplication,
                                      GatherDelly.dellyInsertion,
                                      GatherDelly.dellyInversion,
                                      GatherLumpy.lumpyVCF,
                                      Manta.mantaVCF,                
                                    ]
        
        Array[File] VCFs = select_all(potentialVCF)

        if (runSURVIVOR) {
            call SURVIVOR as vcfSurvivor {
                input:
                    vcfs = VCFs,
                    bamBase = P2Prep.bamBase,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }

        if (runJasmine) {
            call Jasmine as vcfJasmine {
                input:
                    vcfs = VCFs,
                    bamBase = P2Prep.bamBase,
                    runtime_zones = runtime_zones,
                    max_retries = max_retries,
                    preemptible_tries = preemptible_tries
            }
        }
    }

    if ((runSURVIVOR || runJasmine) && visualizeSVs) {
        if (runSURVIVOR) {
            File survivorVCF = select_first([ svtypedSurvivor.survivorQualVCF, vcfSurvivor.survivorQualVCF ])

            scatter (chromosome in chromosomes) {
                call svvizChrom as vizSURVIVORchrom {
                    input:
                        inputBam = P2Prep.outputBam,
                        inputBai = P2Prep.outputBai,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        vcf = survivorVCF,
                        contig = chromosome,
                        runtime_zones = runtime_zones,
                        max_retries = max_retries,
                        preemptible_tries = preemptible_tries
                }
            }

            call GatherSvviz as vizSURVIVOR {
                input:
                    tars = vizSURVIVORchrom.svvizTar,
                    bamBase = P2Prep.bamBase,
                    merger = "survivor",
                    runtime_zones = runtime_zones,
                    preemptible_tries = preemptible_tries
            }
        }

        if (runJasmine) {
            File jasmineVCF = select_first([ svtypedJasmine.jasmineVCF, vcfJasmine.jasmineVCF ])

            scatter (chromosome in chromosomes) {
                call svvizChrom as vizJasmineChrom {
                    input:
                        inputBam = P2Prep.outputBam,
                        inputBai = P2Prep.outputBai,
                        refFasta = refFasta,
                        refIndex = refIndex,
                        vcf = jasmineVCF,
                        contig = chromosome,
                        runtime_zones = runtime_zones,
                        max_retries = max_retries,
                        preemptible_tries = preemptible_tries
                }
            }

            call GatherSvviz as vizJasmine {
                input:
                    tars = vizJasmineChrom.svvizTar,
                    bamBase = P2Prep.bamBase,
                    merger = "jasmine",
                    runtime_zones = runtime_zones,
                    preemptible_tries = preemptible_tries
            }
        }
    }

    output {
        File? breakdancerCTX = GatherBreakdancer.breakdancerCTX
        File? breakdancerVCF = GatherBreakdancer.breakdancerVCF
        File? bdTyped = typeBreakdancer.svtypedVCF

        File? breakseqGFF = Breakseq.breakseqGFF
        File? breakseq_genotypedGFF = Breakseq.breakseq_genotypedGFF
        File? breakseqVCF = Breakseq.breakseqVCF
        File? breakseqVCFindex = Breakseq.breakseqVCFindex
        File? breakseqBAM = Breakseq.breakseqBAM
        File? bsTyped = typeBreakseq.svtypedVCF

        File? CNVnatorOutput = GatherCNVnator.CNVnatorOutput
        File? CNVnatorVCF = GatherCNVnator.CNVnatorVCF
        File? CNVnTyped = typeCNVnator.svtypedVCF

        File? dellyDel = GatherDelly.dellyDeletion
        File? dellyDup = GatherDelly.dellyDuplication
        File? dellyIns = GatherDelly.dellyInsertion
        File? dellyInv = GatherDelly.dellyInversion
        File? dDelTyped = typeDellyDel.svtypedVCF
        File? dDupTyped = typeDellyDup.svtypedVCF
        File? dInsTyped = typeDellyIns.svtypedVCF
        File? dInvTyped = typeDellyInv.svtypedVCF

        File? lumpyGFF = GatherLumpy.lumpyGFF
        File? lumpyVCF = GatherLumpy.lumpyVCF
        File? lumpyTyped = typeLumpy.svtypedVCF

        File? mantaVCF = Manta.mantaVCF
        File? mantaStats = Manta.mantaStats
        File? mantaVariants = Manta.mantaVariants
        File? mantaTyped = typeManta.svtypedVCF

        File? survivorTypedSortedVCF = svtypedSurvivor.survivorSortedVCF
        File? survivorTypedQualVCF = svtypedSurvivor.survivorQualVCF
        File? survivorSortedVCF = vcfSurvivor.survivorSortedVCF
        File? survivorQualVCF = vcfSurvivor.survivorQualVCF

        File? jasmineTyped = svtypedJasmine.jasmineVCF
        File? jasmineUntyped = vcfJasmine.jasmineVCF

        File? survivorViz = vizSURVIVOR.finalSvviz
        File? jasmineViz = vizJasmine.finalSvviz
    }
}

###
# GENERAL TASKS
###
# P2Prep: Generates contigs and bamBase string for later use
task P2Prep {
    input {
        File inputFile
        File inputIndex
        File refFasta
        Boolean filterContigs
        Int preemptible_tries
        String runtime_zones 
    }

    String inputName = "~{basename(inputFile)}"
    String bamName =  sub("~{inputName}", ".cram", ".bam")

    command <<<
        bam_base="~{inputName}"
        if [[ "${bam_base}" == *".cram" ]]; then
            samtools view -b -@ "$(nproc)" -T "~{refFasta}" -o "~{bamName}" "~{inputFile}"
            samtools view -H "~{bamName}" | python /opt/bin/get_contigs.py "~{filterContigs}" > contigs
            samtools index "~{bamName}"
            bam_base="${bam_base%.cram}"
        else
            samtools view -H "~{inputFile}" | python /opt/bin/get_contigs.py "~{filterContigs}" > contigs
            mv "~{inputFile}" "~{bamName}"
            mv "~{inputIndex}" "~{bamName}.bai"
            bam_base="${bam_base%.bam}"
        fi

        ref_genome_type=""
        if [[ "~{basename(refFasta)}" == *"hg19"* ]]; then
            ref_genome_type="hg19"
        elif [[ "~{basename(refFasta)}" == *"hg38"* || "~{basename(refFasta)}" == *"GRCh38"* || "~{basename(refFasta)}" == *"hs38DH"* || "~{basename(refFasta)}" == *"38.fasta" ]]; then
            ref_genome_type="hg38"
        elif [[ "~{basename(refFasta)}" == *"hs37d5"* ]]; then
            ref_genome_type="hs37d5"
        else
            ref_genome_type="other"
        fi

        echo "${ref_genome_type}"
        echo "${bam_base}"
    >>>

    Int diskGb = ceil(10.0 * size(inputFile, "G"))
    
    runtime {
        docker : "szarate/p2_prep:v0.0.1"
        disks : "local-disk ${diskGb} SSD"
        memory: "10G"
        cpu : 8
        zones: runtime_zones
        preemptible: preemptible_tries
    }

    output {
        File contigs = "contigs"
        File outputBam = "~{bamName}"
        File outputBai = "~{bamName}.bai"
        String refType = read_lines(stdout())[0]
        String bamBase = read_lines(stdout())[1]
    }
}

# SplitByChrom: Splits a BAM file into a given chromosome using sambamba
task SplitByChrom {
    input {
        File inputBam
        File inputBai
        String contig
        Int max_retries
        Int preemptible_tries
        String runtime_zones 
    }

    command <<<
        sambamba view -h -f bam -t "$(nproc)" "~{inputBam}" "~{contig}" > "~{contig}.bam"
        sambamba index -t "$(nproc)" "~{contig}.bam"
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "mgibio/sambamba-cwl:0.6.4"
        disks : "local-disk ${diskGb} SSD"
        memory: "10G"
        cpu : 8
        zones: runtime_zones
        maxRetries: max_retries
        preemptible: preemptible_tries
    }

    output {
        File splitBam = "~{contig}.bam"
        File splitBai = "~{contig}.bam.bai"
    }
}

###
# SV CALLERS
###
# Breakdancer
task PrepareBreakdancer {
    input {
        File inputBam
        Int preemptible_tries
        String runtime_zones 
    }

    command <<<
        bam2cfg -x 5000 -o breakdancer.cfg "~{inputBam}"
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "szarate/breakdancer:v1.4.3"
        disks : "local-disk ${diskGb} SSD"
        memory: "5 GB"
        cpu: 1
        zones: runtime_zones
        preemptible: preemptible_tries
    }

    output {
        File configFile = "breakdancer.cfg"
    }
}

task BreakdancerChrom {
    input {
        File inputBam
        File inputBai
        File configFile
        String contig
        Int max_retries
        Int preemptible_tries
        String runtime_zones 
    }

    command <<<
        breakdancer-max "~{configFile}" "~{inputBam}" -o ~{contig} > breakdancer-~{contig}.ctx
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "szarate/breakdancer:v1.4.3"
        disks : "local-disk ${diskGb} SSD"
        memory: "8G"
        cpu : 1
        zones: runtime_zones
        maxRetries: max_retries
        preemptible: preemptible_tries
    }

    output {
        File breakdancerOutput = "breakdancer-~{contig}.ctx"
    }
}

task GatherBreakdancer {
    input {
        Array[File] breakdancerCtx
        String bamBase
        String runtime_zones
        Int preemptible_tries
    }

    command <<<
        cat ~{sep=' ' breakdancerCtx} > "~{bamBase}.breakdancer.ctx"

        python /opt/bin/merge_files.py 1.0 "~{bamBase}.breakdancer.ctx" "~{bamBase}"
        python /opt/bin/ctx_to_vcf.py < "~{bamBase}.breakdancer.ctx" > "~{bamBase}.breakdancer.vcf"
    >>>

    runtime {
        docker : "szarate/breakdancer:v1.4.3"
        memory: "5G"
        cpu : 1
        zones: runtime_zones
        preemptible: preemptible_tries
    }

    output {
        File breakdancerCTX = "~{bamBase}.breakdancer.ctx"
        File breakdancerVCF = "~{bamBase}.breakdancer.vcf"
    }
}

# BreakSeq
task Breakseq {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
        String reference
        String runtime_zones
        Int max_retries
        Int preemptible_tries
    }

    String bamBase='~{basename(inputBam,".bam")}'
    String breakpointLibrary = if (reference == "hg19") then "/breakseq2_bplib_20150129.hg19/breakseq2_bplib_20150129.hg19.gff" else (if reference == "hg38" then "/bplib.hg38.gff" else (if reference == "hs37d5" then "/breakseq2_bplib_20150129.hs37d5/breakseq2_bplib_20150129.gff" else ""))

    command <<<
        mkdir -p "breakseq2"
        #gunzip "~{refFasta}"

        #refName="~{refFasta}"
        #refName="${refName%.gz}"

        /miniconda/bin/run_breakseq2.py \
            --reference "~{refFasta}" \
            --bams "~{inputBam}" \
            --work breakseq2 \
            --bwa /miniconda/bin/bwa \
            --samtools /miniconda/bin/samtools \
            --bplib_gff "~{breakpointLibrary}" \
            --nthreads "$(nproc)" \
            --sample "~{bamBase}"

            gunzip breakseq2/breakseq.vcf.gz

            mv breakseq2/breakseq.vcf "~{bamBase}.breakseq.vcf"
            mv breakseq2/breakseq.vcf.gz.tbi "~{bamBase}.breakseq.vcf.gz.tbi"
            mv breakseq2/breakseq.gff "~{bamBase}.breakseq.gff"
            mv breakseq2/breakseq_genotyped.gff "~{bamBase}_genotyped.breakseq.gff"
            mv breakseq2/final.bam "~{bamBase}.breakseq.bam"
    >>>
    
    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "szarate/breakseq2:v2.2"
        disks : "local-disk ${diskGb} SSD"
        cpu : 8
        zones: runtime_zones
        maxRetries: max_retries
        preemptible: preemptible_tries
    }

    output {
        File breakseqGFF = "${bamBase}.breakseq.gff"
        File breakseq_genotypedGFF = "${bamBase}_genotyped.breakseq.gff"
        File breakseqVCF = "${bamBase}.breakseq.vcf"
        File breakseqVCFindex = "${bamBase}.breakseq.vcf.gz.tbi"
        File breakseqBAM = "${bamBase}.breakseq.bam"
    }
}

# CNVnator
task CNVnatorChrom {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
        String contig
        String runtime_zones
        Int max_retries
        Int preemptible_tries
    }

    command <<<
        cnvnator -root output.root"~{contig}" -chrom "~{contig}" -genome "~{refFasta}" -tree "~{inputBam}"
        cnvnator -root output.root"~{contig}" -chrom "~{contig}" -genome "~{refFasta}" -his 100
        cnvnator -root output.root"~{contig}" -chrom "~{contig}" -genome "~{refFasta}" -stat 100
        cnvnator -root output.root"~{contig}" -chrom "~{contig}" -genome "~{refFasta}" -partition 100
        cnvnator -root output.root"~{contig}" -chrom "~{contig}" -genome "~{refFasta}" -call 100  > "output.cnvnator_calls-~{contig}"
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "szarate/cnvnator:v0.4.1"
        disks : "local-disk ${diskGb} SSD"
        cpu : 8
        zones: runtime_zones
        maxRetries: max_retries
        preemptible: preemptible_tries
    }

    output {
        File CNVnatorOutput = "output.cnvnator_calls-~{contig}"
    }
}

task GatherCNVnator {
    input {
        Array[File] CNVnatorCalls
        String bamBase
        String runtime_zones
        Int preemptible_tries
    }

    command <<<
        cat ~{sep=' ' CNVnatorCalls} > "~{bamBase}.cnvnator.output"

        perl /opt/bin/cnvnator2vcf.pl "~{bamBase}.cnvnator.output" > "~{bamBase}.cnvnator.vcf"
    >>>

    runtime {
        docker : "szarate/cnvnator:v0.4.1"
        memory: "5G"
        cpu : 1
        zones: runtime_zones
        preemptible: preemptible_tries
    }

    output {
        File CNVnatorOutput = "~{bamBase}.cnvnator.output"
        File CNVnatorVCF = "~{bamBase}.cnvnator.vcf"
    }
}

# Delly
task DellyChrom {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
        String runtime_zones
        Int preemptible_tries
        Int max_retries
        Int? cpu_override
        Int cpus = select_first([cpu_override, 8])
    }
    
    String contig = '~{basename(inputBam,".bam")}'

    command <<<
        delly call -t DEL -o "~{contig}.delly.deletion.bcf" -g "~{refFasta}" "~{inputBam}"
        delly call -t DUP -o "~{contig}.delly.duplication.bcf" -g "~{refFasta}" "~{inputBam}"
        delly call -t INS -o "~{contig}.delly.insertion.bcf" -g "~{refFasta}" "~{inputBam}"
        delly call -t INV -o "~{contig}.delly.inversion.bcf" -g "~{refFasta}" "~{inputBam}"

        bcftools view "~{contig}.delly.deletion.bcf" > "~{contig}.delly.deletion.vcf"
        bcftools view "~{contig}.delly.duplication.bcf" > "~{contig}.delly.duplication.vcf"
        bcftools view "~{contig}.delly.insertion.bcf" > "~{contig}.delly.insertion.vcf"
        bcftools view "~{contig}.delly.inversion.bcf" > "~{contig}.delly.inversion.vcf"
    >>>

    runtime {
        docker : "szarate/delly:v0.8.3"
        cpu : cpus
        zones: runtime_zones
        preemptible: preemptible_tries
        maxRetries: max_retries
    }

    output {
        File dellyDeletion = "~{contig}.delly.deletion.vcf"
        File dellyDuplication = "~{contig}.delly.duplication.vcf"
        File dellyInsertion = "~{contig}.delly.insertion.vcf"
        File dellyInversion = "~{contig}.delly.inversion.vcf"
    }
}

task GatherDelly {
    input {
        Array[File] dellyDelVCFs
        Array[File] dellyDupVCFs
        Array[File] dellyInsVCFs
        Array[File] dellyInvVCFs
        String bamBase
        String runtime_zones
        Int preemptible_tries
    }

    command <<<
        python /opt/bin/convert_header.py "~{bamBase}" "~{sep=' ' dellyDelVCFs}" | vcf-sort -c | uniq > "~{bamBase}.delly.deletion.vcf"
        python /opt/bin/convert_header.py "~{bamBase}" "~{sep=' ' dellyDupVCFs}" | vcf-sort -c | uniq > "~{bamBase}.delly.duplication.vcf"
        python /opt/bin/convert_header.py "~{bamBase}" "~{sep=' ' dellyInsVCFs}" | vcf-sort -c | uniq > "~{bamBase}.delly.insertion.vcf"
        python /opt/bin/convert_header.py "~{bamBase}" "~{sep=' ' dellyInvVCFs}" | vcf-sort -c | uniq > "~{bamBase}.delly.inversion.vcf"
    >>>

    runtime {
        docker : "szarate/delly:v0.8.3"
        memory: "5G"
        cpu : 1
        zones: runtime_zones
        preemptible: preemptible_tries
    }

    output {
        File dellyDeletion = "~{bamBase}.delly.deletion.vcf"
        File dellyDuplication = "~{bamBase}.delly.duplication.vcf"
        File dellyInsertion = "~{bamBase}.delly.insertion.vcf"
        File dellyInversion = "~{bamBase}.delly.inversion.vcf"
    }
}

# Lumpy
task LumpyChrom {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
        String runtime_zones
        Int preemptible_tries
    }

    String contig = '~{basename(inputBam,".bam")}'
    String refBase = basename(refFasta)
    String lumpyExclude = if (refBase == "*b37*") then "-x /opt/b37.bed" else (if refBase == "*hg38*" then "-x /opt/hg38.bed" else (if refBase == "*hg19" then "-x /opt/hg19.bed" else ""))
   

    command <<<
        lumpyexpress -B "~{inputBam}" -o "lumpy.~{contig}.vcf" "~{lumpyExclude}" -k
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G") + size(refFasta, "G") + 10) 

    runtime {
        docker : "szarate/lumpy-sv:v0.3.0"
        cpu : 8
        disks : "local-disk ${diskGb} SSD"
        zones: runtime_zones
        preemptible: preemptible_tries
    }

    output {
        File lumpyOutput = "lumpy.~{contig}.vcf"
    }
}

task GatherLumpy {
    input {
        Array[File] lumpyVCFs
        String bamBase
        String runtime_zones
        Int preemptible_tries
    }

    command <<<
        python /opt/bin/convert_header.py "~{bamBase}" "~{sep = ' ' lumpyVCFs}" | vcf-sort -c | uniq > "~{bamBase}.lumpy.vcf"

        python /opt/bin/vcf2bedpe.py -i "~{bamBase}.lumpy.vcf" -o "~{bamBase}.lumpy.gff"
        python /opt/bin/lumpy2merge.py "~{bamBase}.lumpy.gff" "${prefix}" 1.0

        ls -sh
    >>>
    
    
    runtime {
        docker : "szarate/lumpy-sv:v0.3.0"
        memory: "5G"
        cpu : 1
        zones: runtime_zones
        preemptible: preemptible_tries
    }

    output {
        File lumpyGFF = "~{bamBase}.lumpy.gff"
        File lumpyVCF = "~{bamBase}.lumpy.vcf"
    }
}

# Manta
task Manta {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
        File contigs
        String runtime_zones
        Int preemptible_tries
        Int max_retries
    }

    String bamBase='~{basename(inputBam,".bam")}'

    command <<<
        mkdir -p "manta"
        gunzip "~{refFasta}"

        refName="~{refFasta}"
        refName="${refName%.gz}"
        region_string=""

        while read line; do
            region_string="$region_string --region=$line"
        done < "~{contigs}"

        python /usr/local/bin/configManta.py --referenceFasta "${refName}" --normalBam "~{inputBam}" --runDir manta $region_string

        python manta/runWorkflow.py -m local -j "$(nproc)"

        tar -czf stats.tar.gz -C manta/results/stats/ .
        tar -czf variants.tar.gz -C manta/results/variants/ .

        gunzip manta/results/variants/diploidSV.vcf.gz
        mv manta/results/variants/diploidSV.vcf "~{bamBase}.manta.vcf"
        mv stats.tar.gz "~{bamBase}_stats.manta.vcf.gz"
        mv variants.tar.gz "~{bamBase}_variants.manta.vcf.gz"
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "szarate/manta:v1.6.0"
        disks : "local-disk ${diskGb} SSD"
        cpu : 8
        zones: runtime_zones
        preemptible: preemptible_tries
        maxRetries: max_retries
    }

    output {
        File mantaVCF = "~{bamBase}.manta.vcf"
        File mantaStats = "~{bamBase}_stats.manta.vcf.gz"
        File mantaVariants = "~{bamBase}_variants.manta.vcf.gz"
    }
}

### POST-SV CALLING TASKS
# SVTyper: genotypes SVs
task SVTyper {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
        File inputVCF
        String runtime_zones
        Int preemptible_tries
        Int max_retries
        Int? cpu_override
        Int cpus = select_first([cpu_override, 8])
    }

    String vcfBase='~{basename(inputVCF,".vcf")}'

    command <<<
        python /opt/bin/ciend_cpos.py < "~{inputVCF}" 1000 > tmp.vcf
        svtyper-sso --core "$(nproc)" -i tmp.vcf --bam "~{inputBam}" --ref_fasta "~{refFasta}" -o "~{vcfBase}.svtyped.vcf"
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "szarate/svtyper:v0.7.0"
        disks : "local-disk ${diskGb} SSD"
        cpu : cpus
        zones: runtime_zones
        preemptible: preemptible_tries
        maxRetries: max_retries
    }

    output {
        File svtypedVCF = "~{vcfBase}.svtyped.vcf"
    }
}

# SURVIVOR: gathers files
task SURVIVOR {
    input {
        Array[File] vcfs
        String bamBase
        String runtime_zones
        Int preemptible_tries
        Int max_retries
    }
    
    command <<<
        echo "~{sep='\n' vcfs}" > survivor_inputs
        SURVIVOR merge survivor_inputs 1000 1 1 0 0 10 "~{bamBase}.vcf"

        vcf-sort -c < "~{bamBase}.vcf" > "~{bamBase}.survivor.sorted.vcf"

        sed -i 's/SAMPLE/breakdancer/g' "~{bamBase}.survivor.sorted.vcf"
        python3 /opt/bin/combine_combined.py "~{bamBase}.survivor.sorted.vcf" "~{bamBase}" survivor_inputs /opt/bin/all.phred.txt | python3 /opt/bin/correct_max_position.py > "~{bamBase}.survivor.qual.vcf"
    >>>

    runtime {
        docker : "szarate/survivor:1d1d33b"
        memory: "5G"
        cpu : 1
        zones: runtime_zones
        preemptible: preemptible_tries
        maxRetries: max_retries
    }

    output {
        File survivorSortedVCF = "~{bamBase}.survivor.sorted.vcf"
        File survivorQualVCF = "~{bamBase}.survivor.qual.vcf"
    }
}

# Jasmine: alternative to SURVIVOR
task Jasmine {
    input {
        Array[File] vcfs
        String bamBase
        String runtime_zones
        Int preemptible_tries
        Int max_retries
    }

    command <<<
        echo "~{sep='\n' vcfs}" > jasmine_inputs
        java -jar /Jasmine/jasmine.jar -cp src Main file_list=jasmine_inputs out_file="~{bamBase}.jasmine.vcf"
    >>>

    runtime {
        docker : "szarate/jasmine:v1.0.6"
        memory: "5G"
        cpu : 1
        zones: runtime_zones
        preemptible: preemptible_tries
        maxRetries: max_retries
    }

    output {
        File jasmineVCF = "~{bamBase}.jasmine.vcf"
    }
}

# svviz: visualize SVs
task svvizChrom {
    input {
        File inputBam
        File inputBai
        File refFasta
        File refIndex
        File? vcf
        String contig
        String runtime_zones
        Int max_retries
        Int preemptible_tries
    }

    command <<<
        grep "\#" "~{vcf}" > header.txt
        grep "^~{contig}" "~{vcf}" > "~{contig}".txt
        cat header.txt "~{contig}".txt > "~{contig}".vcf

        mkdir -p "~{contig}_outputs/"
        svviz --pair-min-mapq 30 --max-deletion-size 5000 --max-reads 10000 --fast --type batch --summary svviz_summary.tsv -b "~{inputBam}" "~{refFasta}" "~{contig}.vcf" --export "~{contig}_outputs/"

        tar -czvf "~{contig}.tar.gz" "~{contig}_outputs/"
    >>>

    Int diskGb = ceil(2.0 * size(inputBam, "G"))

    runtime {
        docker : "szarate/svviz:v1.6.2"
        disks : "local-disk ${diskGb} SSD"
        cpu : 8
        zones: runtime_zones
        preemptible: preemptible_tries
        maxRetries: max_retries        
    }

    output {
        File svvizTar = "~{contig}.tar.gz"
    }
}

task GatherSvviz {
    input {
        Array[File] tars
        String bamBase
        String merger
        String runtime_zones
        Int preemptible_tries        
    }

    command <<<
        mkdir -p svviz_outputs/
        echo "~{sep='\n' tars}" > all_tars

        while read line; do
            tar -xzf "${line}" -C svviz_outputs/
        done < all_tars
        tar -czf "~{bamBase}.~{merger}.svviz.tar.gz" svviz_outputs/
    >>>
    
    runtime {
        docker : "ubuntu:20.04"
        memory: "5G"
        cpu : 1
        zones: runtime_zones
        preemptible: preemptible_tries
    }

    output {
        File finalSvviz = "~{bamBase}.~{merger}.svviz.tar.gz"
    }
}