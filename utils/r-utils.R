## list here all the main variables and functions used my scripts
## if any of these are used more than once in different scripts then they will be sourced from here
library(yaml) 
library(RColorBrewer)
library(data.table)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(dplyr)

set.seed(2022)
project_config <- suppressWarnings(read_yaml('./config/project_config.yaml')) ## this because genome size is generally > max integer range printed by R. leave it like this unless u'll need this entry

##=============##
## DIRECTORIES ##
##=============##
base_dir = project_config$project_config$baseDir 
project_name = project_config$project_config$project 
project_dir = paste(base_dir,'/',project_name,'/',sep='')
vast_dir = project_config$project_config$tmp_dir
data_dir = paste(project_dir,'data/',sep='')
outdir = project_config$outdir
qcdir = project_config$qcdir
tables_dir = paste(outdir,'tables/',sep='')
plots_dir = paste(outdir,'plots/',sep='')
bam_dir = paste(outdir, project_config$align,sep='') 
logs_dir = project_config$logs
peaks_dir = paste(outdir,project_config$`peak-calling`,sep='')
consensus_peaks_dir = paste(outdir,project_config$`consensus-peak`,sep='')

## setwd for all R scripts
setwd(project_dir)

##===========##
## VARIABLES ##
##===========##
species = project_config$species$species
standard_chr = unlist(lapply(project_config$standard_chromosomes,function(x)paste0("chr",x,sep='')))
grange_cols = c('seqnames','start','end') 

## if your species is not any of these two then change these lines accordingly
if (species %like% 'musculus'){
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  transcripts <- GenomicFeatures::transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
} else if (species %like% 'sapiens'){
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  transcripts <- GenomicFeatures::transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
}

samples = project_config$samples

qualitative_palette = brewer.pal.info[brewer.pal.info$category == 'qual',]
sample_palette = sample(unlist(mapply(brewer.pal, qualitative_palette$maxcolors, rownames(qualitative_palette))),length(samples))
names(sample_palette) = samples

##===========##
## FUNCTIONS ##
##===========##
## function that creates a dir with specified path if it doesnt exist (atm is just like dir.create - will change)
create_dir <- function(path){
  dir <- path
  dir.create(dir,showWarnings=F,recursive = T)
  return(dir)
}
