#!/usr/bin/env python3

import multiprocessing
import pandas
from pathlib import Path


def demux_target(wildcards):
    '''
    glob the demux output directory to see which samples worked
    '''
    cdir = checkpoints.demultiplex.get(**wildcards).output[0]
    my_samples = glob_wildcards(Path(cdir, '{sample}_r1.fastq.gz')).sample
    output_dict = {
        'files': expand('output/000_tmp/reads/{sample}_r{r}.fastq.gz',
                        sample=my_samples,
                        r=['1', '2']),
        'directory': cdir}
    return(output_dict)


def resolve_path(x):
    return(Path(x).resolve().as_posix())


###########
# GLOBALS #
###########

hyp_ref = 'data/mhyp.fa'
aeth_ref = 'data/maeth.fa'
sample_csv = 'data/samples.csv'

# software
bbmap = 'shub://TomHarrop/seq-utils:bbmap_38.76'
csdemux = 'shub://TomHarrop/csdemux:csdemux_v0.0.4'
honeybee_genotype_pipeline = (
    'shub://TomHarrop/'
    'honeybee-genotype-pipeline:honeybee_genotype_pipeline_v0.0.11')
plink = 'shub://MarissaLL/singularity-containers:plink_1.9'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'


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
               ref=['hyp', 'aeth']),
        expand('output/020_filtered/{ref}/pruned.stats.txt',
               ref=['hyp', 'aeth'])

# fst stats etc.
# plink --vcf pruned.vcf.gz --fst --double-id --allow-extra-chr --within within.txt
# plink --vcf pruned.vcf.gz --het --double-id --allow-extra-chr --within within.txt
# paste \
#     <( bcftools query -l pruned.vcf.gz | cut -d'_' -f2 )
#     <( bcftools query -l pruned.vcf.gz) \
#     <( bcftools query -l pruned.vcf.gz | cut -d'_' -f2,3 )
# paste blah2 blah2 blah1 > within.txt

# prune LD with bcftools
rule prune_vcf:
    input:
        vcf = 'output/020_filtered/{ref}/filtered.vcf.gz'
    output:
        temp('output/020_filtered/{ref}/pruned.vcf')
    log:
        'output/logs/prune_vcf.{ref}.log'
    singularity:
        samtools
    shell:
        'bcftools +prune '
        '--max-LD 0.1 '
        '--window 50 '
        '{input.vcf} '
        '> {output} '
        '2> {log}'


rule filter_vcf:
    input:
        vcf = 'output/010_genotypes/{ref}/calls.vcf.gz',
    output:
        temp('output/020_filtered/{ref}/filtered.vcf')
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


rule genotype:
    input:
        csv = 'output/005_config/samples.csv',
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

# try manual demux
rule generate_sample_csv:
    input:
        unpack(demux_target),
        sample_csv = 'data/samples.csv'
    output:
        csv = 'output/005_config/samples.csv'
    script:
        'src/generate_sample_csv.py'

checkpoint demultiplex:
    input:
        r1 = 'data/muxed/Undetermined_S0_L006_R1_001.fastq.gz',
        r2 = 'data/muxed/Undetermined_S0_L006_R2_001.fastq.gz',
        sample_csv = 'data/samples.csv'
    output:
        directory('output/010_demux/reads')
    log:
        'output/logs/demultiplex.log'
    params:
        outdir = 'output/010_demux'
    threads:
        workflow.cores
    singularity:
        csdemux
    shell:
        'csdemux '
        '--threads {threads} '
        '--mem_gb 500 '
        '--restart_times 2 '
        '--r1 {input.r1} '
        '--r2 {input.r2} '
        '--samples_csv {input.sample_csv} '
        '--outdir {params.outdir} '
        '&> {log}'


# generic vcf index
rule generic_index_vcf:
    input:
        Path('{folder}', '{file}.vcf')
    wildcard_constraints:
        folder = 'output/(?!010).*'
    output:
        gz = Path('{folder}', '{file}.vcf.gz'),
        tbi = Path('{folder}', '{file}.vcf.gz.tbi')
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} '
        '; '
        'tabix -p vcf {output.gz}'


# generic vcf stats
rule generic_vcf_stats:
    input:
        vcf = Path('{folder}', '{file}.vcf.gz'),
        tbi = Path('{folder}', '{file}.vcf.gz.tbi')
    wildcard_constraints:
        folder = 'output/(?!010).*'
    output:
        stats = Path('{folder}', '{file}.stats.txt'),
        plot = directory(Path('{folder}', '{file}.plots'))
    singularity:
        samtools
    shell:
        'bcftools stats '
        '--verbose '
        '-S <( bcftools query -l {input.vcf} ) '
        '{input.vcf} '
        '> {output.stats} ; '
        'plot-vcfstats -s -v -p {output.plot} {output.stats} '
