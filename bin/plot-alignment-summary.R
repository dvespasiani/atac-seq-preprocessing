## use this script to plot a summary of Bowtie2 alignment results
library(dplyr)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

source('./utils/r-utils.R')

encode_acceptable = 0.8
encode_ideal = 0.95

input_file = unlist(snakemake@input[[1]])
output_plot = unlist(snakemake@output[[1]])
pattern = paste(paste(samples,'.log',sep=''),collapse='|')
files <- list.files(dirname(input_file),recursive=F,full.names=T,pattern=pattern)

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

qc_results <- Map(mutate,qc_results,sample=samples)%>%rbindlist()

cat('Plotting the Bowtie2 alignment summary results \n')

plot_align <- function(results,ylabel){
    p <- ggplot(results,aes(x=sample,y=fraction_reads,fill=sample))+ 
        geom_bar(stat="identity")+xlab('')+ylab(ylabel)+
        ylim(0,1)+
        scale_fill_manual(values=sample_palette)+
        theme_classic()+
        theme(
            legend.position='none',
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
            )
    return(p)
}

plot_alignment_rate <- plot_align(
    qc_results[class %in% 'alignment_rate'],'Bowtie2 alignment rate')+ 
    geom_hline(yintercept=encode_acceptable)+
    geom_hline(yintercept=encode_ideal,linetype="dashed")+
    geom_text(aes(1.5,encode_acceptable,label = 'acceptable', vjust = -1))+
    geom_text(aes(1.5,encode_ideal,label = 'ideal', vjust = -1))

plot_mappable_rate <- plot_align(qc_results[class %in% 'uniquely_mappable'],'Uniquely mappable rate')

pdf(output_plot,width=7,height=7)
grid.arrange(plot_alignment_rate, plot_mappable_rate, nrow =2)
dev.off()

