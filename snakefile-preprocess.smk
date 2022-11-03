import os
import itertools
import glob

##===================##
##   config params   ##
##===================##
configfile: "config/snakemake-config.yaml"

species = config['species']
assembly = config['assembly']

basedir = config['basedir']
fastqdir = config['fastqdir']

adapters = config['adapters']
read_length = config['read_length']
genome_size = config['genome_size']
genome_index = config['genome_index']
blacklist = config['blacklist']
chrom_sizes = config['chrom_sizes']

sample = config['samples']

# genome2bit_index = config['genome2bit_index']
# merged_sample = config['merged_sample']

# sample = list(itertools.chain(*config['augmented_samples'])) ## flattens list of lists

##=================##
## I/O directories ##
##=================##
outdir = 'out/preprocessing/'
qcdir = outdir + 'qc/'
logs = 'logs/'

##=================##
##   Parameters    ##
##=================##
read_minQ = config['read_minQ']
npeaks = config['npeaks']
fragment_size = config['fragment_size']
shift = config['shift']
pval_thresh = config['pval_thresh']

##=================##
##   rule groups   ##
##=================##
main = 'main'
qc = 'qc'

##=================##
##  include rules  ##
##=================##

include: "rules/fastqc.smk"
include: "rules/trim-adapter.smk"
include: "rules/alignment.smk"
include: "rules/post-alignment.smk"
include: "rules/peak-calling.smk"
include: "rules/deeptools.smk"


rule all:
      input:  
            ## fastqc
            expand(fastqdir + "{sample}_R1_001.fastq.gz",sample=sample),
            expand(fastqdir + "{sample}_R2_001.fastq.gz",sample=sample),
            expand(outdir + "fastqc/{sample}_R1_001_fastqc.zip", sample=sample),
            expand(outdir + "fastqc/{sample}_R2_001_fastqc.zip", sample=sample),
            
            ## trim-adapter
            expand(outdir + "trim-adapter/{sample}-1-trimmed.fastq.gz",sample=sample),
            expand(outdir + "trim-adapter/{sample}-2-trimmed.fastq.gz",sample=sample),
            expand(outdir + "trim-adapter/{sample}-1-unpaired.fastq.gz", sample=sample),
            expand(outdir + "trim-adapter/{sample}-2-unpaired.fastq.gz", sample=sample),
            
            ## alignment
            expand(outdir + "alignment/{sample}.bam",sample=sample),

            ## post-alignment
            expand(outdir + "post-alignment/{sample}-nochrM.bam",sample=sample),
            expand(outdir + "post-alignment/{sample}-nochrM-encodefiltered.bam",sample=sample),
            expand(outdir + "post-alignment/{sample}-nochrM-encodefiltered-fixmate.bam",sample=sample),
            expand(outdir + "post-alignment/{sample}-nochrM-encodefiltered-fixmate-rmorphanread.bam",sample=sample),
            expand(outdir + "post-alignment/{sample}-nochrM-encodefiltered-fixmate-rmorphanread-dupmark.bam",sample=sample),
            expand(qcdir + "{sample}-duplicate-rate.qc",sample=sample),
            expand(outdir + "post-alignment/{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bam",sample=sample),
            expand(outdir + "post-alignment/{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bai",sample=sample),
            expand(outdir + "post-alignment/{sample}-tn5-shifted.bam",sample=sample),
            expand(outdir + "post-alignment/{sample}-tn5-shifted-sorted.bam",sample = sample),
            expand(outdir + "post-alignment/{sample}-tn5-shifted-sorted.bam.bai",sample = sample),

            ## peak calling
            expand(outdir + "peak_calling/{sample}-macs2_peaks.xls",sample=sample),
            expand(outdir + "peak_calling/{sample}-macs2_treat_pileup.bdg",sample=sample),
            expand(outdir + "peak_calling/{sample}-macs2_summits.bed",sample=sample),
            expand(outdir + "peak_calling/{sample}-macs2_control_lambda.bdg",sample=sample),
            expand(outdir + "peak_calling/{sample}-macs2_peaks.narrowPeak",sample=sample),
            expand(outdir + "peak_calling/{sample}-macs2-peaks-filtered.narrowPeak.gz",sample=sample),
            expand(outdir + "peak_calling/{sample}-macs2-peaks-filtered-sorted.narrowPeak.gz",sample=sample),
            # expand(outdir + "peak_calling/{sample}-macs2_FE.bdg",sample=sample),
            # expand(outdir + "peak_calling/{sample}-fe-signal.bedgraph",sample=sample),
            # expand(outdir + "peak_calling/{sample}-fe-signal-sorted.bedgraph",sample=sample),
            # expand(outdir + "peak_calling/{sample}-fe-signal.bigwig",sample=sample),
            # expand(outdir + "peak_calling/{sample}-ppois-sval",sample=sample),
            # expand(outdir + "peak_calling/{sample}-macs2_ppois.bdg",sample=sample),
            # expand(outdir + "peak_calling/{sample}-ppois-signal.bedgraph",sample=sample),
            # expand(outdir + "peak_calling/{sample}-ppois-signal-sorted.bedgraph",sample=sample),
            # expand(outdir + "peak_calling/{sample}-ppois-signal.bigwig",sample=sample),
            expand(qcdir + "{sample}-frip.txt",sample=sample),

            
            ## deeptools
            expand(outdir + "deeptools/{sample}-noblacklist.bam",sample=sample),
            expand(outdir + "deeptools/{sample}-noblacklist.bai",sample=sample),
            expand(outdir + "deeptools/{sample}-SeqDepthNorm.bw",sample=sample),
            outdir  + "deeptools/samples-bam-coverage.png",
            outdir + "deeptools/samples-plot-fingerprint.png",
            outdir + "deeptools/multiBAM-fingerprint-metrics.txt",
            outdir + "deeptools/multiBAM-fingerprint-rawcounts.txt",
            expand(outdir + "deeptools/{sample}-GC-content.txt",sample=sample),
            expand(outdir + "deeptools/{sample}-plot-GC-content.png",sample=sample),
            outdir + "deeptools/multibam-summary.npz",
            outdir + "deeptools/multibam-readcounts.txt",
            outdir + "deeptools/pearson-corr-multibam.png",
            outdir + "deeptools/pearson-corr-multibamsum-matrix.txt"
