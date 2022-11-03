#!/usr/bin/env Rscript

## use this script to plot a summary of Bowtie2 alignment results
library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)

source('./utils/r-utils.R')

encode_acceptable = 0.2
encode_ideal = 0.3

outplot_dir <- create_dir(path=paste(plots_dir,'atac-seq-qc',sep='')) 
input_dir = paste(logs_dir,'alignment',sep='')

files <- list.files(input_dir,full.names=T,recursive=F)

cat("Parsing all these files:", paste(basename(files),collapse=' , '), " located in the ", input_dir, " directory \n")

qc_results <- lapply(files,function(x) {
    x <- fread(x,sep='\t',col.names='bowtie_output')[
        c(3,14),
        ][
            ,oar:=ifelse(bowtie_output %like% 'overall',
            gsub("\\%.*","",bowtie_output), 
            gsub("\\%.*","",gsub("[\\(\\)]", "", regmatches(bowtie_output, gregexpr("\\(.*?\\)", bowtie_output)))))
            ][
                ,fraction_reads:=as.numeric(oar)/100
                ][
                    ,class:=ifelse(bowtie_output %like% 'overall','alignment_rate','uniquely_mappable')
                    ][
                        ,c('fraction_reads','class')
                        ]
})

sample_names <- gsub(".*/","",gsub('\\.log.*','',files))
qc_results <- Map(mutate,qc_results,sample=sample_names)%>%rbindlist()

cat('Plotting the Bowtie2 alignment summary results into the ', outplot_dir, ' directory \n')

pdf(paste(outplot_dir,'bowtie2-alignment-qc.pdf',sep=''),width=7,height=7)
ggplot(qc_results,aes(x=sample,y=fraction_reads,fill=sample))+ 
geom_bar(stat="identity")+xlab('')+ylab('Bowtie2 Alignment rate')+
geom_hline(yintercept=encode_acceptable)+
scale_fill_manual(values=sample_palette)+
geom_text(aes(1.5,encode_acceptable,label = 'acceptable', vjust = -1))+
geom_text(aes(1.5,encode_ideal,label = 'ideal', vjust = -1))+
geom_hline(yintercept=encode_ideal,linetype="dashed")+
theme_classic()+
theme(
    legend.position='none',
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
dev.off()

cat('Done \n')
