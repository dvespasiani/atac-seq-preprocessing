import sys
sys.path.append('/home/users/allstaff/vespasiani.d/.local/lib/python3.6/site-packages')
import yaml
import subprocess

config_file = './config/snakemake-config.yaml'

with open(config_file, 'r') as file:
    config = yaml.safe_load(file)  

binpath = "./bin/"

def run_executable(exe, i = None,o = None):
    command = binpath + exe
    if i is None:
        return (subprocess.call ([command]))
    elif o is None:
        return (subprocess.call ([command, "-i " + i ]))
    else:
        return (subprocess.call ([command, "-i " + i , "-o " + o]))

## count number reads in each bam file 
run_executable(exe = 'get-counts.sh', i = config['outdir'] + 'post-alignment' , o = config['qcdir']+ "qc-number-reads-bams.txt")

## count number peaks in each narrowPeak file 
run_executable(exe = 'get-counts.sh', i = config['outdir'] + 'peak_calling' , o = config['qcdir']+ "qc-number-peaks-narrowPeaks.txt")

## library complexity
run_executable(exe = 'estimate-lib-complexity.sh', i = config['outdir'] + 'post-alignment' , o = config['qcdir'])

## alignment summary
run_executable(exe = 'plot-alignment-summary.R')

## frip
run_executable(exe = 'plot-frip-summary.R')

## peakqc
run_executable(exe = 'plot-peak-qcs.R')

## tss enrich
run_executable(exe = 'plot-tss-enrich.R')
