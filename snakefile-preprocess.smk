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
# chrom_sizes = config['chrom_sizes']

sample = config['samples']


# genome2bit_index = config['genome2bit_index']
# merged_sample = config['merged_sample']

# augmented_samples = list(itertools.chain(*config['augmented_samples'])) ## flattens list of lists

##=====================
## I/O directories
##=====================
outdir = 'out/preprocessing/'
logs = 'logs/'

# aligment_dir =  {'output': 'output/Alignment/Files/','logs':'output/Alignment/logs/','qc':'output/Alignment/qc/'}
# aligment_outdir = aligment_dir['output']
# aligment_logdir = aligment_dir['logs']
# aligment_qcdir = aligment_dir['qc']

# peakcall_dir = {'output': 'output/PeakCalling/Files/','logs':'output/PeakCalling/logs/','qc':'output/PeakCalling/qc/'}
# peakcall_outdir = peakcall_dir['output']
# peakcall_logdir = peakcall_dir['logs']
# peakcall_qcdir = peakcall_dir['qc']

# deeptools_dir = {'output': 'output/Deeptools/Files/','logs':'output/Deeptools/logs/','qc':'output/Deeptools/qc/'}
# deeptools_outdir = deeptools_dir['output']
# deepTools_logdir = deeptools_dir['logs']
# deeptools_qcdir = deeptools_dir['qc']

##===========
## params
##===========
read_minQual = 20

##================
## rule groups
##================
main = 'main'
qc = 'qc'
# deeptools_group = 'deeptools'

##================
## include rules
##================
include: "rules/fastqc.smk"
include: "rules/trim-adapter.smk"
include: "rules/alignment.smk"
include: "rules/post-alignment.smk"
# include: "rules/peakcall.smk"
# include: "rules/cross_corr.smk"
# include: "rules/deeptools.smk"


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
            expand(outdir + "post-alignment/{sample}-dup.qc",sample=sample),
            expand(outdir + "post-alignment/{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bam",sample=sample),
            expand(outdir + "post-alignment/{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bai",sample=sample),
            expand(outdir + "post-alignment/{sample}-tn5-shifted.bam",sample=sample),
            expand(outdir + "post-alignment/{sample}-tn5-shifted-sorted.bam",sample = sample),
            expand(outdir + "post-alignment/{sample}-tn5-shifted-sorted.bam.bai",sample = sample),

            
## Peak call
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_peaks.xls",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_summits.bed",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_treat_pileup.bdg",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_control_lambda.bdg",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_default_peaks.narrowPeak",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_macs2_default_peaks_filtered.narrowPeak",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_macs2_default_peaks_filtered_sorted.narrowPeak.gz",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_macs2_FE.bdg",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_fc_signal.bedgraph",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_fc_signal_sorted.bedgraph",augmented_samples=augmented_samples),
            expand(peakcall_outdir + "{augmented_samples}_fc_signal.bigwig",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_ppois_sval",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_macs2_ppois.bdg",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_ppois_signal.bedgraph",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_ppois_signal_sorted.bedgraph",augmented_samples=augmented_samples),
            # expand(peakcall_outdir + "{augmented_samples}_ppois_signal.bigwig",augmented_samples=augmented_samples),
            expand(peakcall_qcdir + "{augmented_samples}_default.frip.txt",augmented_samples=augmented_samples),
            # "output/PeakCalling/ConsensusPeaks/consensus_peak.bed.gz"

## Cross correlation
      # expand("output/ENCODE_CC/Files/{sample}_1_trimmed.fastq.gz",sample=sample),
      # expand("output/ENCODE_CC/Files/{sample}_1_trimmed.bam",sample=sample),
      # expand("output/ENCODE_CC/Files/{sample}_1_trimmed_q30.bam",sample=sample),
      # expand("output/ENCODE_CC/Files/{sample}_R1_trimmed_q30_SE.tagAlign.gz",sample=sample),
      # expand("output/ENCODE_CC/Files/{sample}.filt.sample.25Mreads.SE.tagAlign.gz",sample=sample),
      # expand("output/ENCODE_CC/QCs/{sample}_filt_25Mreads.SE.cc.qc",sample=sample),
      # expand("output/ENCODE_CC/QCs/{sample}_filt_25Mreads.SE.cc.plot.pdf",sample=sample),

## Deeptools 
            # expand(deeptools_outdir + "{sample}.SeqDepthNorm.bw",sample=sample),
            # deeptools_outdir + "Samples_plotCoverage.png",
            # deeptools_outdir + "Samples_plotfingerprint.png",
            # deeptools_outdir + "multiBAM_fingerprint_metrics.txt",
            # deeptools_outdir + "multiBAM_fingerprint_rawcounts.txt",
            # expand(deeptools_outdir + "{sample}_GC_content.txt",sample=sample),
            # expand(deeptools_outdir + "{sample}_plot_GC_content.png",sample=sample),
            # deeptools_outdir + "Summary.npz",
            # deeptools_outdir + "Readcounts.txt",
            # deeptools_outdir + "PearsonCor_multibamsum.png",
            # deeptools_outdir + "PearsonCor_multibamsum_matrix.txt"
