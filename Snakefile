#!/usr/bin/env python3

import multiprocessing
import pandas
from pathlib import Path


def resolve_path(x):
    return(Path(x).resolve().as_posix())


###########
# GLOBALS #
###########

hyp_ref = 'data/mhyp.fa'
aeth_ref = 'data/maeth.fa'
sample_csv = 'data/samples.csv'

# software
honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.10')
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'
plink = 'shub://MarissaLL/singularity-containers:plink_1.9'

########
# MAIN #
########

# get a list of samples
sample_data = pandas.read_csv(sample_csv,
                              index_col='sample')
all_samples = sorted(set(sample_data.index))

#########
# RULES #
#########

rule target:
    input:
        # expand('output/010_genotypes/{ref}/calls.vcf.gz',
        #        ref=['hyp', 'aeth'])
        expand('output/020_filtered/{ref}/pruned.vcf.gz',
               ref=['hyp', 'aeth'])

# get a set of LD-free SNPs
rule prune_vcf:
    input:
        vcf = 'output/020_filtered/{ref}/filtered.vcf',
        prune = 'output/020_filtered/{ref}/plink.prune.in'
    output:
        vcf = 'output/020_filtered/{ref}/pruned.vcf'
    log:
        'output/logs/prune_vcf.{ref}.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '-i \'ID=@{input.prune}\' '
        '{input.vcf} '
        '> {output.vcf} '
        '2> {log}'

rule list_pruned_snps:
    input:
        vcf = 'output/020_filtered/{ref}/filtered.vcf'
    output:
        'output/020_filtered/{ref}/plink.prune.in'
    params:
        vcf = lambda wildcards, input: resolve_path(input.vcf),
        wd = 'output/020_filtered/{ref}',
        indep = '50 10 0.1'     # 50 kb window, 10 SNPs, r2 < 0.1
    log:
        resolve_path('output/logs/list_pruned_snps.{ref}.log')
    singularity:
        plink
    shell:
        'cd {params.wd} || exit 1 ; '
        'plink '
        '--vcf {params.vcf} '
        '--double-id '
        '--allow-extra-chr '
        '--set-missing-var-ids @:# '
        '--indep-pairwise {params.indep} '
        '--out plink '
        '&> {log}'


rule filter_vcf:
    input:
        vcf = 'output/010_genotypes/{ref}/calls.vcf.gz',
    output:
        'output/020_filtered/{ref}/filtered.vcf'
    params:
        min_maf = 0.05,
        f_missing = 0.2
    log:
        'output/logs/filter_vcf.{ref}.log'
    singularity:
        samtools
    shell:
        'bcftools view '
        '--min-af {params.min_maf}:nonmajor '
        '--exclude "F_MISSING>{params.f_missing}" '
        # '-S <( cut -f1 {input.popmap} ) '
        '{input.vcf} '
        '> {output} '
        '2> {log}'


checkpoint genotype:
    input:
        csv = sample_csv,
        ref = lambda wildcards: hyp_ref if wildcards.ref == 'hyp' else aeth_ref
    output:
        # disable until the pipeline works again
        # bam = 'output/010_genotypes/{ref}/merged.bam',
        # cutoffs = 'output/010_genotypes/{ref}/040_stats/ldepth.mean_cutoffs.csv',
        # fai = 'output/010_genotypes/{ref}/015_ref/ref.fasta.fai',
        # ref = 'output/010_genotypes/{ref}/015_ref/ref.fasta',
        vcf = 'output/010_genotypes/{ref}/calls.vcf.gz',
    params:
        wd = 'output/010_genotypes/{ref}',
        ploidy = '2'
    log:
        'output/logs/genotype.{ref}.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        honeybee_genotype_pipeline
    shell:
        'honeybee_genotype_pipeline '
        '--ref {input.ref} '
        '--samples_csv {input.csv} '
        '--outdir {params.wd} '
        '--ploidy {params.ploidy} '
        '--threads {threads} '
        '--restart_times 1 '
        '&> {log}'



# generic vcf index
rule generic_index_vcf:
    input:
        Path('{folder}', '{file}.vcf')
    wildcard_constraints:
        file = '(?!calls).*'
    output:
        gz = Path('{folder}', '{file}.vcf.gz'),
        tbi = Path('{folder}', '{file}.vcf.gz.tbi')
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} '
        '; '
        'tabix -p vcf {output.gz}'
