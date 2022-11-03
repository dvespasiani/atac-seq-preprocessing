#!/usr/bin/env bash
## original script comes from the ENCODE https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit
## input: filtered bams (dup removed)
## output: table with different metrics used to estimate the library complexity

helpFunction()
{
   echo ""
   echo "Usage: $0 -i input_dir -o out_dir"
   echo -e "\t-i Directory where the bams/narrowPeak files are stored"
   echo -e "\t-o path for the output directory"
   exit 1 # Exit script after printing help
}

while getopts "i:o:" flag; do
    case "${flag}" in
        i) 
            input_dir="$OPTARG"
            ;;
        o) 
            out_dir="$OPTARG"
            ;;
        ?) 
            helpFunction 
            ;; 
    esac
done
shift $((OPTIND -1))

tmp_file="${out_dir}tmp_file.qc"
tmp_filename="${out_dir}tmp_filename.txt"
final_outfile="${out_dir}library-complexity.txt"

echo 
echo
echo "Estimating library complexity for filtered bam files located in the $input_dir directory"
echo
echo


for f in ${input_dir}/*-nodup.bam ; do

    echo
    echo "Estimating library complexity for $f"
    echo
    echo "$(echo $(basename $f)| cut -f 1 -d '-')" > ${tmp_filename}
    bedtools bamtobed -i ${f} | \
    awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | \
    grep -v 'chrM' | sort | uniq -c | \
    awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} \
    END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > results.txt

    paste -d "\t" ${tmp_filename} results.txt >> ${tmp_file}
done 

echo -e "sample\tTotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF\tPBC1\tPBC2" | cat - ${tmp_file} > ${final_outfile}

rm ${tmp_file} && rm results.txt && rm ${tmp_filename}

echo
echo
echo "Done, the output file is this one here $final_outfile"
echo
echo

