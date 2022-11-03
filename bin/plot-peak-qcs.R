#!/usr/bin/env Rscript

## Script to produce some extra QCs on MACS2 peaks 
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)

source('./utils/r-utils.R')

outplot_dir <- create_dir(path=paste(plots_dir,'atac-seq-qc',sep='')) 

pattern = "*-filtered-sorted.narrowPeak.gz"
input_dir = paste(outdir,'peak_calling/',sep='')

cat ('Parsing all the ', pattern, 'files, located in the ' , input_dir, '\n')

peak_files <- list.files(input_dir,full.names=T, recursive=F,pattern=pattern)

peaks <- lapply(peak_files,function(x)
    x<-fread(x,sep='\t',header=F,select=c(1:3),col.names=grange_cols)[
            ,width:=end-start
        ]
)
names(peaks) = samples
peaks <- Map(mutate,peaks,samples=samples)%>%rbindlist()

peaks <- peaks[,rounded_width:=plyr::round_any(width, 100)] 

cat('Plotting the peak QC results into the ', outplot_dir, ' directory \n')

pdf(paste(outplot_dir,'distribution-peak-sizes.pdf',sep=''),width=7,height=10)
ggplot(peaks,aes(x=rounded_width,fill=samples))+
    geom_bar()+
    scale_fill_manual(values=sample_palette)+
    facet_wrap(samples~.,ncol=2)+
    ylab('Number of peaks')+ xlab('Peak size')+
    theme_classic()+
    theme(
        legend.position='none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
        )
dev.off()

cat('Done \n')
