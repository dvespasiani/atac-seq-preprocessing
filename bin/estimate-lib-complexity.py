#!/usr/bin/env python

## original script comes from the ENCODE https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
## use this script to create a table with different metrics used to estimate the library complexity
## input files are all filtered and duplicates removed bam files

import os
from collections import defaultdict
import pandas as pd
import subprocess

# input_file = "out/alignment/all_samples-tn5-shifted-sorted.bam"
input_file = snakemake.input[0]
output = snakemake.output[0]

input_dir = os.path.dirname(input_file) + '/'
files = [input_dir + f for f in os.listdir(input_dir) if f.endswith('fixmate.bam')]

keys = ['sample','TotalReadPairs','DistinctReadPairs','OneReadPair','TwoReadPairs','NRF','PBC1','PBC2']
list_dict = []

for f in files:
    fname = os.path.splitext(os.path.basename(f))[0]
    fname = fname.split('-nochrM', 1)[0]
    cmd = "bedtools bamtobed -i " + f + """| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c |""" + """ awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t""%d\t""%d\t""%d\t""%f\t""%f\t""%f\t",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' """
    libcomplx = subprocess.check_output(cmd, shell=True)
    libcomplx = libcomplx.decode("utf-8")
    libcomplx = libcomplx.split()
    libcomplx.insert(0, fname)
    list_dict.append(dict(zip(keys,libcomplx)))

res = defaultdict(list)

for d in list_dict:
    for k, v in d.items():
        res[k].append(v)

res = dict(res)
df = pd.DataFrame.from_dict(res)

df.to_csv(output,index=False,header=True,sep='\t')

print('')
print('Done')
print('')

