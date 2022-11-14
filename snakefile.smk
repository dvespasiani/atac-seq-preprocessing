import os
import itertools
import glob

##===================##
##   config params   ##
##===================##
configfile: "config/project_config.yaml"

##=================##
## I/O directories ##
##=================##
project_config = config['project_config']
project_name = project_config['project']
inputdir = project_config['datadir']
outdir = project_config['outdir']
qcdir = project_config['qcdir']
logs = project_config['logs']
plots = project_config['plots']
tables = project_config['tables']

##===========##
##  Species  ##
##===========##
species_config = config['species_info']
species = species_config['species']
assembly = species_config['assembly']
genome_size = species_config['genome_size']
genome_index = species_config['genome_index']
blacklist = species_config['blacklist']
genome2bit_index = species_config['genome2bit_index']
chrom_sizes = species_config['chrom_sizes']

##=====================##
##   ATAC experiment   ##
##=====================##
atac_specifics = config['atac_specifics']
adapters = atac_specifics['adapters']
read_length = atac_specifics['read_length']

##==============##
##   Samples    ##
##==============##
sample = config['samples']
all_samples = config['all_samples']
combined_sample = list(itertools.chain(*config['combined_sample'])) ## this flattens a list of lists

##================##
##   QC params    ##
##================##
read_minQ = config['read_minQ']
npeaks = config['npeaks']
fragment_size = config['fragment_size']
shift = config['shift']
pval_thresh = config['pval_thresh']

##=================##
##     Rules       ##
##=================##
main = 'main'
qc = 'qc'

rulename_fastqc = config['fastqc'] +'/'
rulename_trim = config['trim-adapter'] +'/'
rulename_alignment = config['align'] +'/'
rulename_peak = config['peak-calling'] +'/'
rulename_deeptools = config['deeptools'] +'/'
rulename_conspeak = config['consensus-peak'] + '/'
rulename_qc = config['atac-qc'] + '/'

include: "rules/fastqc.smk"
include: "rules/trim-adapter.smk"
include: "rules/alignment.smk"
include: "rules/peak-calling.smk"
include: "rules/deeptools.smk"
include: "rules/qcs.smk"
include: "rules/consensus-peak.smk"

rule all:
      input:  
            ## fastqc
            expand(inputdir + "{sample}_R1_001.fastq.gz",sample=sample),
            expand(inputdir + "{sample}_R2_001.fastq.gz",sample=sample),
            expand(outdir + rulename_fastqc + "{sample}_R1_001_fastqc.zip", sample=sample),
            expand(outdir + rulename_fastqc + "{sample}_R2_001_fastqc.zip", sample=sample),
            
            ## trim-adapter
            expand(outdir  + rulename_trim + "{sample}-1-trimmed.fastq.gz",sample=sample),
            expand(outdir  + rulename_trim + "{sample}-2-trimmed.fastq.gz",sample=sample),
            expand(outdir  + rulename_trim + "{sample}-1-unpaired.fastq.gz", sample=sample),
            expand(outdir  + rulename_trim + "{sample}-2-unpaired.fastq.gz", sample=sample),
            
            # ## alignment
            expand(outdir  + rulename_alignment + "{sample}.bam",sample=sample),
            expand(outdir  + rulename_alignment + "{sample}-nochrM.bam",sample=sample),
            expand(outdir  + rulename_alignment + "{sample}-nochrM-encodefiltered.bam",sample=sample),
            expand(outdir  + rulename_alignment + "{sample}-nochrM-encodefiltered-fixmate.bam",sample=sample),
            expand(outdir  + rulename_alignment + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread.bam",sample=sample),
            expand(outdir  + rulename_alignment + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread-dupmark.bam",sample=sample),
            expand(qcdir + rulename_alignment + "{sample}-duplicate-rate.qc",sample=sample),
            expand(outdir  + rulename_alignment + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bam",sample=sample),
            expand(outdir  + rulename_alignment + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bai",sample=sample),
            expand(outdir  + rulename_alignment + "{sample}-tn5-shifted.bam",sample=sample),
            expand(outdir  + rulename_alignment + "{all_samples}-tn5-shifted.bam",all_samples = all_samples),
            expand(outdir  + rulename_alignment + "{combined_sample}-tn5-shifted-sorted.bam",combined_sample = combined_sample),
            expand(outdir  + rulename_alignment + "{combined_sample}-tn5-shifted-sorted.bam.bai",combined_sample = combined_sample),

            ## peak calling
            expand(outdir + rulename_peak + "{combined_sample}-macs2_peaks.xls",combined_sample=combined_sample),
            expand(outdir + rulename_peak + "{combined_sample}-macs2_treat_pileup.bdg",combined_sample=combined_sample),
            expand(outdir + rulename_peak + "{combined_sample}-macs2_summits.bed",combined_sample=combined_sample),
            expand(outdir + rulename_peak + "{combined_sample}-macs2_control_lambda.bdg",combined_sample=combined_sample),
            expand(outdir + rulename_peak + "{combined_sample}-macs2_peaks.narrowPeak",combined_sample=combined_sample),
            expand(outdir + rulename_peak + "{combined_sample}-macs2-peaks-filtered.narrowPeak.gz",combined_sample=combined_sample),
            expand(outdir + rulename_peak + "{combined_sample}-macs2-peaks-filtered-sorted.narrowPeak.gz",combined_sample=combined_sample),
            # expand(outdir + rulename_peak + "{combined_sample}-macs2_FE.bdg",combined_sample=combined_sample),
            # expand(outdir + rulename_peak + "{combined_sample}-fe-signal.bedgraph",combined_sample=combined_sample),
            # expand(outdir + rulename_peak + "{combined_sample}-fe-signal-sorted.bedgraph",combined_sample=combined_sample),
            # expand(outdir + rulename_peak + "{combined_sample}-fe-signal.bigwig",combined_sample=combined_sample),
            # expand(outdir + rulename_peak + "{combined_sample}-ppois-sval",combined_sample=combined_sample),
            # expand(outdir + rulename_peak + "{combined_sample}-macs2_ppois.bdg",combined_sample=combined_sample),
            # expand(outdir + rulename_peak + "{combined_sample}-ppois-signal.bedgraph",combined_sample=combined_sample),
            # expand(outdir + rulename_peak + "{combined_sample}-ppois-signal-sorted.bedgraph",combined_sample=combined_sample),
            # expand(outdir + rulename_peak + "{combined_sample}-ppois-signal.bigwig",combined_sample=combined_sample),
            expand(qcdir + "{combined_sample}-frip.txt",combined_sample=combined_sample),

            
            ## deeptools
            expand(outdir + rulename_deeptools + "{sample}-noblacklist.bam",sample=sample),
            expand(outdir + rulename_deeptools + "{sample}-noblacklist.bai",sample=sample),
            expand(outdir + rulename_deeptools + "{sample}-SeqDepthNorm.bw",sample=sample),
            outdir  + rulename_deeptools + "samples-bam-coverage.png",
            outdir + rulename_deeptools + "samples-plot-fingerprint.png",
            outdir + rulename_deeptools + "multiBAM-fingerprint-metrics.txt",
            outdir + rulename_deeptools + "multiBAM-fingerprint-rawcounts.txt",
            expand(outdir + rulename_deeptools + "{sample}-GC-content.txt",sample=sample),
            expand(outdir + rulename_deeptools + "{sample}-plot-GC-content.png",sample=sample),
            outdir + rulename_deeptools + "multibam-summary.npz",
            outdir + rulename_deeptools + "multibam-readcounts.txt",
            outdir + rulename_deeptools + "pearson-corr-multibam.png",
            outdir + rulename_deeptools + "pearson-corr-multibamsum-matrix.txt",

            # ## qcs
            tables + rulename_qc + 'number-reads-bam-files.txt',
            tables + rulename_qc + 'number-peaks.txt',
            tables + rulename_qc + "library-complexity.txt",
            plots + rulename_qc + "bowtie2-alignment-summary.pdf",
            plots + rulename_qc + 'frip-summary-qc.pdf',
            plots + rulename_qc + 'extra-peak-qcs.pdf',
            plots + rulename_qc + 'tss-enrichment.pdf',

            # consensus peak
            plots + rulename_qc + 'support-consensus-peak.pdf',
            outdir + rulename_conspeak + 'consensus-peak.bed.gz',

