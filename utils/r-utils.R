## list here all the main variables and functions used my scripts
## if any of these are used more than once in different scripts then they will be sourced from here

##=============##
## DIRECTORIES ##
##=============##
base_dir = '/stornext/General/data/user_managed/grpu_jchoi_0/projects/davide/'
project = 'atac-pipeline'
vast_dir = paste('/vast/scratch/users/vespasiani.d/',project,'/',sep='')
project_dir = paste(base_dir,project,'/',sep='')
data_dir = paste(project_dir,'data/',sep='')
out_dir = paste(project_dir,'out/',sep='')
tables_dir = paste(out_dir,'tables/',sep='')
plots_dir = paste(out_dir,'plots/',sep='')

##===========##
## VARIABLES ##
##===========##
species = 'musmusculus'
standard_chr <- paste0("chr", c(1:20),sep='')

##===========##
## FUNCTIONS ##
##===========##
## function that creates a dir with specified path if it doesnt exist (atm is just like dir.create - will change)
create_dir <- function(path){
  dir <- path
  dir.create(dir,showWarnings=F,recursive = T)
  return(dir)
}
