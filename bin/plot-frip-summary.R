library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(argparse)

setwd('/stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/atac-pipeline')

source('./utils/r-utils.R')

encode_acceptable = 0.8
encode_ideal = 0.95

outplot_dir <- create_dir(path=paste(plots_dir,'atac-seq-qc',sep='')) 

parser <- ArgumentParser()
parser$add_argument("-i", "--input_dir",help="Input directory containing the Bowtie2 log files")
args <- parser$parse_args()

cat('Parsing all files in input directory \n')
## FRiP
files <- list.files(args$i,full.names=T,recursive=F,pattern='frip.txt')
qc_results <- lapply(files,function(x){
    x <- fread(x,sep='\t',header=T,col.names='results')[2,][
        ,frip:=round(as.numeric(gsub('.*=','',results)),2)
        ]
})

sample_names <- gsub(".*/","",gsub('\\-frip.txt*','',files))
qc_results <- Map(mutate,qc_results,sample=sample_names)%>%rbindlist()
qc_results <- qc_results[,c('frip','sample')][,type:='reads_in_peak']


cat('Plotting FRiP summary results \n')

pdf(paste(outplot_dir,'frip-qc.pdf',sep=''),width=7,height=7)
ggplot(qc_results,aes(x=sample,y=frip,fill=sample))+ 
geom_bar(stat='identity')+xlab('')+ylab('FRiP score')+
geom_hline(yintercept=encode_acceptable)+
geom_text(aes(1.5,encode_acceptable,label = 'acceptable', vjust = -1))+
geom_hline(yintercept=encode_ideal,linetype="dashed")+
geom_text(aes(1.5,encode_ideal,label = 'ideal', vjust = -1))+
# scale_fill_manual(values=samples_palette)+
theme_classic()+
theme(
    legend.position='none',
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
    )
dev.off()

cat('Done \n')
