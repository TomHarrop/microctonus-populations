#!/usr/bin/env python3

import pandas
from pathlib import Path
import snakemake as sm

print(snakemake.input)

sample_data = pandas.read_csv(snakemake.input['sample_csv'],
                              index_col='sample')

cdir = snakemake.input['directory']
my_bcs = sm.io.glob_wildcards(Path(cdir, '{bc}_r1.fastq.gz')).bc
my_mask = sample_data['barcode'].isin(my_bcs)
my_csv = sample_data[my_mask]
my_csv['r1_path'] = [f'output/000_tmp/reads/{x}_r1.fastq.gz'
                     for x in my_csv['barcode']]
my_csv['r2_path'] = [f'output/000_tmp/reads/{x}_r2.fastq.gz'
                     for x in my_csv['barcode']]
my_csv.to_csv(snakemake.output['csv'])
