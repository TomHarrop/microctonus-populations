#!/usr/bin/env python3

import multiprocessing
import pandas

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
samtools = 'shub://TomHarrop/singularity-containers:samtools_1.9'

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
        expand('output/010_genotypes/{ref}/calls.vcf.gz',
               ref=['hyp', 'aeth'])

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