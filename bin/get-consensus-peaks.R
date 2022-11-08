
## Script used to define the set of consensus peaks identified across multiple samples
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(UpSetR)
library(GenomicRanges)

source('./utils/r-utils.R')

sample_peak_files = as.list(snakemake@input[["sample_peaks"]])
combined_peak_file = unlist(snakemake@input[["combined_peaks"]])
output_plot = snakemake@output[[1]]
outfile = snakemake@output[[2]]

sample_peaks <- lapply(sample_peak_files,function(x){
    x <- fread(x,sep='\t',header=F,select=c(1:4),col.names=c(grange_cols,'peakID'))[!seqnames %like% '_']%>%setorderv(grange_cols)%>%unique()
    setkeyv(x,grange_cols)
    return(x)
})

names(sample_peaks) = samples
sample_peaks <- Map(mutate,sample_peaks,sample=samples)

combined_peak <- fread(
    combined_peak_file,sep='\t',header=F,select=c(1:4),col.names=c(grange_cols,'peakID')
)[!seqnames %like% '_'][,sample := snakemake_config$all_samples]%>%setorderv(grange_cols)%>%unique()

setkeyv(combined_peak,grange_cols)

cat("Intersecting the set of individual peaks with the combined file to identify consensus peaks\n")

overlapping_peakIDs <- lapply(
    sample_peaks,function(x){ 
        overlap <- foverlaps(x,combined_peak,type='any')%>%na.omit()
        # return the list of peakIDs of the combined sample that have an overlapping peak in the individual peak sample
        peakIDs <- overlap$peakID
        return(peakIDs)
})

individual_peak_support <- unlist(copy(overlapping_peakIDs))

individual_peak_support <- data.table(peakID = individual_peak_support)
individual_peak_support <- individual_peak_support[,support:=.N,by=.(peakID)]%>%unique()

combined_peak_support <- merge(copy(combined_peak), copy(individual_peak_support), by = "peakID", all = TRUE)
combined_peak_support[is.na(combined_peak_support)] <- 0


cat("Plotting an Upset plot to visualise how many samples support the peaks identified in the combined file \n")

pdf(output_plot,width=7,height=7)
upset(fromList(overlapping_peakIDs),nsets = length(samples), order.by = "freq")
dev.off()

cat("Exporting a consensus peak bed file with the number of samples supporting the peaks identified in the combined file \n")

## export new set of peaks
fwrite(combined_peak_support,outfile,sep='\t',col.names=T,quote=F,row.names=F)
cat("Done \n")
