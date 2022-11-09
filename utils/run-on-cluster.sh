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
modules="./config/modules.txt"
eval "$(cat "$modules")"


if [ ! -d 'out/' ]; then
  mkdir -p 'out/';
fi

if [ ! -d 'logs/' ]; then
  mkdir -p 'logs/';
fi


snakemake --profile ./config/snakemake-slurm-config.yaml 