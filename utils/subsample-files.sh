#!/bin/bash
## This script extract a random number (n) of reads from all the fastq files located in the dir. 
## Subsampled files can then be used to test the Snakemake pipeline 

usage(){
   echo ""
   echo "Usage: $0 -i <indir> -o <outdir> -n <numbreads>"
   echo -e "\t-i the input dir with all the original fastq files"
   echo -e "\t-o the output dir where to save all the subsampled fastq files"
   echo -e "\t-n the number of reads to subsample"
   echo ""
   exit 1 # Exit script after printing help
}

while getopts "i:o:n:" opt; do
    case "${opt}" in
        i) 
            indir=${OPTARG}
            ;;
        o) 
            outdir=${OPTARG}
            ;;
        n) 
            numbreads=${OPTARG}
            ;;
        *)
            echo "invalid command $OPTARG"
            ;;
        ?) 
            usage 
            ;; # Print helpFunction
    esac
done
shift $((OPTIND -1))


## create outdir if it doesnt exists
if [ ! -d $outdir ]; then
  echo
  echo "creating $outdir"
  echo
  mkdir -p $outdir;
fi

## subsample files
for f in $indir/*.fastq*; do
  fname="$(echo $(basename $f)| cut -f 1 -d '/')"
  echo "sampling $numbreads reads from $fname and saving it to $outdir/${fname}"
  echo
  if [[ $f =~ \.gz$ ]]; then
    gzip -cd "$f" | head -n $numbreads  | gzip -nc > "$outdir/${fname}" 
  else
    head $f -n $numbreads  | gzip -nc > "$outdir/${fname}.gz"
  fi
done

