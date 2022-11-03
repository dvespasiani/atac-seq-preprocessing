#!/usr/bin/env Rscript

library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
# library(argparse)

source('./utils/r-utils.R')

encode_acceptable = 0.8
encode_ideal = 0.95

outplot_dir <- create_dir(path=paste(plots_dir,'atac-seq-qc',sep='')) 

## FRiP
files <- list.files(qcdir,full.names=T,recursive=F,pattern='frip.txt')

cat("Parsing all these files:", paste(basename(files),collapse=' , '), " located in the ", qcdir, " directory \n")

qc_results <- lapply(files,function(x){
    x <- fread(x,sep='\t',header=T,col.names='results')[2,][
        ,frip:=round(as.numeric(gsub('.*=','',results)),2)
        ]
})

sample_names <- gsub(".*/","",gsub('\\-frip.txt*','',files))
qc_results <- Map(mutate,qc_results,sample=sample_names)%>%rbindlist()
qc_results <- qc_results[,c('frip','sample')][,type:='reads_in_peak']


cat('Plotting the FRiP summary results into the ', outplot_dir, ' directory \n')

pdf(paste(outplot_dir,'frip-qc.pdf',sep=''),width=7,height=7)
ggplot(qc_results,aes(x=sample,y=frip,fill=sample))+ 
geom_bar(stat='identity')+xlab('')+ylab('FRiP score')+
geom_hline(yintercept=encode_acceptable)+
geom_text(aes(1.5,encode_acceptable,label = 'acceptable', vjust = -1))+
geom_hline(yintercept=encode_ideal,linetype="dashed")+
geom_text(aes(1.5,encode_ideal,label = 'ideal', vjust = -1))+
scale_fill_manual(values=sample_palette)+
theme_classic()+
theme(
    legend.position='none',
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
dev.off()

cat('Done \n')
