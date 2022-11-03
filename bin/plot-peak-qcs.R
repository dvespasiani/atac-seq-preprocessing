## Script to produce some extra QCs on MACS2 peaks 
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)
library(ggpubr)

source('./utils/r-utils.R')

outplot_dir <- create_dir(path=paste(plots_dir,'atac-seq-qc',sep='')) 

cat ('Reading in all filtered and sorted narrowPeaks \n')

peak_files <- list.files(paste(out_dir,'peak_calling/',sep=''),full.names=T, recursive=F,pattern="*-filtered-sorted.narrowPeak.gz")

peaks <- lapply(peak_files,function(x)
    x<-fread(x,sep='\t',header=F,select=c(1:3),col.names=grange_cols)[
            ,width:=end-start
        ]
)
names(peaks) = samples
peaks <- Map(mutate,peaks,samples=samples)%>%rbindlist()

peaks <- peaks[,rounded_width:=plyr::round_any(width, 100)] 


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
