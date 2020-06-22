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
bbmap = 'shub://TomHarrop/seq-utils:bbmap_38.76'

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


def demux_target(wildcards):
    cdir = checkpoints.demultiplex.get(**wildcards).output[0]
    my_bcs = glob_wildcards(Path(cdir, '{bc}_r1.fastq.gz')).bc
    output_dict = {
        'files': expand('output/000_tmp/reads/{barcode}_r{r}.fastq.gz',
                        barcode=my_bcs,
                        r=['1', '2']),
        'directory': cdir}
    return(output_dict)

rule generate_sample_csv:
    input:
        unpack(demux_target),
        sample_csv = 'data/samples.csv'
    output:
        csv = 'output/005_config/samples.csv'
    script:
        'src/generate_sample_csv.py'


# try manual demux
checkpoint demultiplex:
    input:
        r1 = 'data/muxed/Undetermined_S0_L006_R1_001.fastq.gz',
        r2 = 'data/muxed/Undetermined_S0_L006_R2_001.fastq.gz'
    output:
        directory('output/000_tmp/reads')
    log:
        'output/logs/demultiplex.log'
    params:
        names = ','.join(sorted(set(sample_data['barcode']))),
        out = 'output/000_tmp/reads/%_r1.fastq.gz',
        out2 = 'output/000_tmp/reads/%_r2.fastq.gz'
    singularity:
        bbmap
    shell:
        'demuxbyname.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={params.out} '
        'out2={params.out2} '
        'names={params.names} '
        'prefixmode=f '
        '-Xmx100g '
        '2> {log}'


