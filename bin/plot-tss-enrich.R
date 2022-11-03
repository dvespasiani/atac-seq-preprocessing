#!/usr/bin/env Rscript

## script used to calculate the TSS enrichment. Use tn5-shifted-sorted.bam files
library(dplyr)
library(data.table)
library(magrittr)
library(GenomicRanges)
library(ggplot2)
library(ggpubr)

source('./utils/r-utils.R')

encode_ideal = 7
encode_acceptable = 5

outplot_dir <- create_dir(path=paste(plots_dir,'atac-seq-qc',sep='')) 

pattern = "*tn5-shifted-sorted.bam$"
cat('Reading all', pattern, 'bam files located in the', bam_dir, ' directory\n')

## read aligned tn5-shifted bams
param <- csaw::readParam(pe = "both",restrict=standard_chr,max.frag=1000)
bams <-  list.files(bam_dir, recursive = T,full.names = T,pattern=pattern)
alignment <- lapply(bams, function(x)GenomicAlignments::readGAlignments(x))

txs <- transcripts

cat('Calculating the TSS enrichment score \n')

get_tsse <- function(genomic_alignment,transcripts){
    tsse <- ATACseqQC::TSSEscore(genomic_alignment, transcripts)
    range_bp <- 100*(-9:10-.5)
    tsse_df <- data.table(range=range_bp,tsse_enrichment=tsse$values,tsse_score=tsse$TSSEscore)
    return(tsse_df)
}

tss_enrichment <- lapply(alignment,function(x)get_tsse(x,txs))
tss_enrichment <- Map(mutate,tss_enrichment,sample=samples)%>%rbindlist()

cat('Plotting the TSS enrichment score into the ', outplot_dir, ' directory \n')

pdf(paste(outplot_dir,'tss-enrichment.pdf',sep=''),width=7,height=7)
ggplot(tss_enrichment,aes(x=range,y=tsse_enrichment,col=sample))+
    geom_line(size=1)+
    geom_hline(yintercept=encode_acceptable)+
    geom_text(aes(-900,encode_acceptable,label = 'acceptable', vjust = -1),col='black')+
    geom_hline(yintercept=encode_ideal,linetype="dashed")+
    geom_text(aes(-900,encode_ideal,label = 'ideal', vjust = -1),col='black')+
    ylab('Enrichment')+xlab('Distance from TSS')+
    scale_color_manual(values=sample_palette)+
    theme_classic()+
    theme(
        legend.position='bottom',
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)
        )
dev.off()

cat('Done \n')



