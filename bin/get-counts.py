## use this script to create a txt file with: 1) filename and 2) number of peaks/reads

import os
import pandas as pd
from snakemake.io import *
import subprocess

input_file = snakemake.input[0]
output = snakemake.output[0]

input_dir = os.path.dirname(input_file) + '/'
files = [input_dir + f for f in os.listdir(input_dir)]

fname = []
counts = []
for f in files:
    if f.endswith('.bam'):
        cmd = "samtools view " + f + '| wc -l ' 
        fname.append(os.path.splitext(os.path.basename(f))[0])
        count = int(subprocess.check_output(cmd, shell=True))
        counts.append(count)
    elif f.endswith('filtered.narrowPeak.gz'):
        fname.append(os.path.splitext(os.path.basename(f))[0])
        peak_df = pd.read_csv(f,header=None,sep='\t',usecols=[0,3],compression='gzip')
        entries = peak_df[peak_df.columns[0]].count()
        counts.append(entries)

data = {'file': fname, 'counts': counts}
df = pd.DataFrame(data)
df.to_csv(output,index=False,header=True,sep='\t')

print('')
print('Done')
print('')