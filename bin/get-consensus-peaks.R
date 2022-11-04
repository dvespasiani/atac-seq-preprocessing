
## Script used to define the set of consensus peaks identified across multiple samples
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(UpSetR)
library(GenomicRanges)

source('./utils/r-utils.R')

args <- commandArgs(trailingOnly=TRUE) 

input_dir = args[[1]]
output_plot = args[[2]]
outfile = args[[3]]

filename = '-macs2-peaks-filtered-sorted.narrowPeak.gz'
pattern = paste(paste(samples,filename,sep=''),collapse='|')

peak_files <- list.files(input_dir,recursive=F,full.names=T,pattern = pattern)

cat("Reading the individual and combined peak files located in the ", peaks_dir, " directory \n")

sample_peaks <- lapply(peak_files,function(x){
    x <- fread(x,sep='\t',header=F,select=c(1:4),col.names=c(grange_cols,'peakID'))[!seqnames %like% '_']%>%setorderv(grange_cols)%>%unique()
    setkeyv(x,grange_cols)
})
names(sample_peaks) = samples
sample_peaks <- Map(mutate,sample_peaks,sample=samples)

snakemake_config$all_samples

combined_peak <- fread(
    paste(peaks_dir,'/',snakemake_config$all_samples,filename, sep=''),
    sep='\t',header=F,select=c(1:4),col.names=c(grange_cols,'peakID')
)[!seqnames %like% '_'][,sample := snakemake_config$all_samples]%>%setorderv(grange_cols)%>%unique()

setkeyv(combined_peak,grange_cols)

cat("Intersecting the set of individual peaks: ", paste(basename(peak_files),collapse=' , '), ", with the combined file:", paste(snakemake_config$all_samples,filename, sep=''), " to identify consensus peaks\n")

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
cat("Done, the plot can be found in the: ", output_plot,  "while here is the new file: ", outfile, "\n")
