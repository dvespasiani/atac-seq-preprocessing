#!/bin/bash 
#SBATCH --job-name="atac-pipeline"
#SBATCH --time 2-00:00:00
#SBATCH --ntasks=10
#SBATCH --partition=regular
#SBATCH --mem=2000
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vespasiani.d@wehi.edu.au
#SBATCH --output=slurm-report/main-%j.out
#SBATCH --error=slurm-report/main-%j.err

## check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

## Run the job from the directory where it was launched (default)

##The modules to load:
module load R/4.2.1 
module load pandoc/2.3.1 
module load deeptools/3.5.1
module load snakemake/7.12.0 
module load fastqc
module load trimmomatic/0.36
module load bowtie2/2.4.4
module load samtools/1.9
module load picard-tools/2.26.11
module load python/3.7.0
module load macs2/2.2.7.1
module load ucsc-tools/331
module load bedtools/2.26.0

## set up and activate a virtual py env in your HPCScratch area
. mypyenv/bin/activate

if [ ! -d 'out/' ]; then
  mkdir -p 'out/';
fi

if [ ! -d 'logs/' ]; then
  mkdir -p 'logs/';
fi


snakemake --profile ./config/slurm/