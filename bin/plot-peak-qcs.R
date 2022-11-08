## Script to produce some extra QCs on MACS2 peaks 
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)

source('./utils/r-utils.R')

input_file = unlist(snakemake@input[[1]])
output_plot = unlist(snakemake@output[[1]])

pattern = paste(paste(samples,"-macs2-peaks-filtered-sorted.narrowPeak.gz",sep=''),collapse="|")
peak_files <- list.files(dirname(input_file),recursive=F,full.names=T,pattern=pattern)

peaks <- lapply(peak_files,function(x)
    x<-fread(x,sep='\t',header=F,select=c(1:3),col.names=grange_cols)[
            ,width:=end-start
        ]
)
names(peaks)  = samples
peaks <- Map(mutate,peaks,sample=samples)%>%rbindlist()

peaks <- peaks[,rounded_width:=plyr::round_any(width, 100)] 

cat('Plotting the peak QC results \n')

pdf(output_plot,width=7,height=7)
ggplot(peaks,aes(x=rounded_width,fill=sample))+
    geom_bar()+
    facet_wrap(~ sample,ncol=floor(length(samples)/2))+
    scale_fill_manual(values=sample_palette)+
    ylab('Number of peaks')+ xlab('Peak size')+
    theme_classic()+
    theme(
        legend.position='none',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
        )
dev.off()

cat('Done \n')
