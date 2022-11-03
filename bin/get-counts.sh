#!/usr/bin/env bash

## use this script to create a spreadsheet with: 1) filename and 2) number of peaks/reads

helpFunction()
{
   echo ""
   echo "Usage: $0 -i input_dir -o out_file"
   echo -e "\t-i Directory where the bams/narrowPeak files are stored"
   echo -e "\t-o path + name of output file"
   exit 1 # Exit script after printing help
}

while getopts "i:o:" flag; do
    case "${flag}" in
        i) 
            input_dir="$OPTARG"
            ;;
        o) 
            outfile="$OPTARG"
            ;;
        ?) 
            helpFunction 
            ;; 
    esac
done
shift $((OPTIND -1))

echo
echo
echo "Getting filenames and counts for the files located in the $input_dir "
echo
echo

for f in ${input_dir}/* ; do

    if [[ "$f" == *.bam ]] || [[ "$f" == *narrowPeak.gz ]] ; then
        fname="$(echo $(basename $f | cut -f 1 -d '.') )"
        printf '%s\n' "$fname"  >> names.out
    fi
        if [[ "$f" == *.bam ]] ; then
        samtools view "$f" | wc -l >> counts.out 
        
        elif [[ "$f" == *narrowPeak.gz ]]; then
            gunzip -nc "$f" | wc -l >> counts.out
        fi
done

echo
echo
echo "Combining all tmp files and creating the final output file"
echo
echo

paste -d "\t" names.out counts.out > tmp.txt
echo -e "file\tcounts" | cat - tmp.txt > $outfile

rm names.out && rm counts.out  && rm tmp.txt 

echo 
echo "Done, the output file is this one here $outfile"
echo
echo