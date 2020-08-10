#!/usr/bin/env python3

import pandas
from pathlib import Path
import snakemake as sm

sample_data = pandas.read_csv(snakemake.input['sample_csv'],
                              index_col='sample')

cdir = snakemake.input['directory']
my_samples = sm.io.glob_wildcards(Path(cdir, '{sample}_r1.fastq.gz')).sample
my_mask = sample_data['sample'].isin(my_samples)
my_csv = sample_data[my_mask]
my_csv['r1_path'] = [f'output/010_demux/reads/{x}_r1.fastq.gz'
                     for x in my_csv['barcode']]
my_csv['r2_path'] = [f'output/010_demux/reads/{x}_r2.fastq.gz'
                     for x in my_csv['barcode']]
my_csv.to_csv(snakemake.output['csv'])
