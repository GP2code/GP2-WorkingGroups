#!/usr/bin/env python
import glob
from pathlib import Path
from natsort import natsorted, ns

parentdir = Path(srcdir("")).parents[0]
configfile: "configs/config_GRCh38_v1.6_noanno.yml"

rule prepare:
    input: 'CADD_input/{chr}/{shard}.vcf'
    output: temp('CADD_tmp/{chr}/{shard}.prepared.vcf')
    conda: '../envs/CADD.yml'
    threads: 1
    shell:
        '''
        cat {input} \
        | python {config[CADD_workdir]}/src/scripts/VCF2vepVCF.py \
        | sort -k1,1 -k2,2n -k4,4 -k5,5 \
        | uniq > {output}
        '''        

rule prescore:
    input: rules.prepare.output
    output:
        novel=temp('CADD_tmp/{chr}/{shard}.novel.vcf'),
        prescored=temp('CADD_tmp/{chr}/{shard}.pre.tsv')
    conda: '../envs/CADD.yml'
    threads: 1
    shell:
        '''
        # Prescoring
        echo '## Prescored variant file' > {output.prescored};
        if [ -d {config[PrescoredFolder]} ]
        then
            for PRESCORED in $(ls {config[PrescoredFolder]}/*.tsv.gz)
            do
                cat {input} \
                | python {config[CADD_workdir]}/src/scripts/extract_scored.py --header \
                    -p $PRESCORED --found_out={output.prescored}.tmp \
                > {input}.tmp;
                cat {output.prescored}.tmp >> {output.prescored}
                mv {input}.tmp {input};
            done;
            rm {output.prescored}.tmp
        fi
        mv {input} {output.novel}
        '''

rule annotation:
    input: 'CADD_tmp/{chr}/{shard}.novel.vcf'
    output: temp('CADD_tmp/{chr}/{shard}.anno.tsv.gz')
    conda: '../envs/CADD.yml'
    threads: 1
    shell:
        '''
        cat {input} \
        | vep --quiet --cache --offline --dir {config[VEPpath]} \
            --buffer 1000 --no_stats --species homo_sapiens \
            --db_version={config[EnsemblDB]} --assembly {config[GenomeBuild]} \
            --format vcf --regulatory --sift b --polyphen b --per_gene --ccds --domains \
            --numbers --canonical --total_length --vcf --force_overwrite --output_file STDOUT \
        | python {config[CADD_workdir]}/src/scripts/annotateVEPvcf.py \
            -c {config[ReferenceConfig]} \
        | gzip -c > {output}
        '''

rule imputation:
    input: 'CADD_tmp/{chr}/{shard}.anno.tsv.gz'
    output: temp('CADD_tmp/{chr}/{shard}.csv.gz')
    conda: '../envs/CADD.yml'
    threads: 1
    shell:
        '''
        zcat {input} \
        | python {config[CADD_workdir]}/src/scripts/trackTransformation.py -b \
            -c {config[ImputeConfig]} -o {output} --noheader;
        '''

rule score:
    input:
        impute='CADD_tmp/{chr}/{shard}.csv.gz',
        anno='CADD_tmp/{chr}/{shard}.anno.tsv.gz'
    output: temp('CADD_tmp/{chr}/{shard}.novel.tsv')
    conda: '../envs/CADD.yml'
    threads: 1
    shell:
        '''
        python {config[CADD_workdir]}/src/scripts/predictSKmodel.py \
            -i {input.impute} -m {config[Model]} -a {input.anno} \
        | python {config[CADD_workdir]}/src/scripts/max_line_hierarchy.py --all \
        | python {config[CADD_workdir]}/src/scripts/appendPHREDscore.py \
            -t {config[ConversionTable]} > {output};
    
        if [ "{config[Annotation]}" = 'False' ]
        then
            cat {output} | cut -f {config[Columns]} | uniq > {output}.tmp
            mv {output}.tmp {output}
        fi
        '''

rule join:
    input:
        pre='CADD_tmp/{chr}/{shard}.pre.tsv',
        novel='CADD_tmp/{chr}/{shard}.novel.tsv'
    output: 'CADD/{chr}/{shard}.tsv.gz'
    conda: '../envs/CADD.yml'
    threads: 1
    shell:
        '''
        (
        echo "{config[Header]}";
        head -n 1 {input.novel};
        cat {input.pre} {input.novel} \
        | grep -v "^#" \
        | sort -k1,1 -k2,2n -k3,3 -k4,4 || true;
        ) | bgzip -c > {output} 

        tabix -s1 -b2 -e2 --force {output}
        '''

rule CADD_per_chr:
    input:
        lambda wildcards: expand('CADD/{{chr}}/{shard}.tsv.gz',chr= chr_dict.keys(), shard=chr_dict.get(wildcards.chr))
    output: 'CADD/{chr}.tsv.gz'
    conda: '../envs/CADD.yml'
    threads: 1
    shell:
        '''
        (
        echo "{config[Header]}";
        zcat {input} | awk 'NR > 1 && /^#/ {{ next }} 1' \
        | grep -v "^#" \
        | sort -k1,1 -k2,2n -k3,3 -k4,4 || true;
        ) | bgzip -c > {output}

        tabix -s1 -b2 -e2 --force {output}
        '''





