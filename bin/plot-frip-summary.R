library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)

source('./utils/r-utils.R')

encode_acceptable = 0.8
encode_ideal = 0.95

args <- commandArgs(trailingOnly=TRUE) 

input_dir = args[[1]]
output_plot = args[[2]]

## FRiP
files <- list.files(input_dir,full.names=T,recursive=F,pattern='frip.txt')

cat("Parsing all these files:", paste(basename(files),collapse=' , '), " located in the ", input_dir, " directory \n")

qc_results <- lapply(files,function(x){
    x <- fread(x,sep='\t',header=T,col.names='results')[2,][
        ,frip:=round(as.numeric(gsub('.*=','',results)),2)
        ]
})

sample_names <- gsub(".*/","",gsub('\\-frip.txt*','',files))
qc_results <- Map(mutate,qc_results,sample=sample_names)%>%rbindlist()
qc_results <- qc_results[,c('frip','sample')][,type:='reads_in_peak']


cat('Plotting the FRiP summary results \n')

pdf(output_plot,width=7,height=7)
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

cat('Done, plot can be found at:' , output_plot, '\n')
