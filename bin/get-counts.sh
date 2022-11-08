input=$(realpath $1)
output=$(realpath $2)

echo "input file" "$input"
echo
input_dir=$(echo "${input%/*}/")
echo "input dir" "$input_dir"
echo

fname="$(echo $(basename $input | cut -f 1 -d '.') )"
"$fname"  > "$output"

# for f in "$input_dir"* ; do
#     fname="$(echo $(basename $f | cut -f 1 -d '.') )"
#     printf '%s\n' "$fname"  > $output
# done
    # if [[ "$f" == *.bam ]] || [[ "$f" == *narrowPeak.gz ]] ; then
    #     fname="$(echo $(basename $f | cut -f 1 -d '.') )"
    #     printf '%s\n' "$fname"  >> names.out
    # fi
        if [[ "$f" == *.bam ]] ; then
        counts=$(samtools view "$f" | wc -l)
        
        elif [[ "$f" == *narrowPeak.gz ]]; then
            counts=$(gunzip -nc "$f" | wc -l)
        fi
    entry=$(paste -d "\t" $fname $counts)
done

echo "entries"
echo $entry
echo

printf '%s\n' "$entry"  > $output



# rm names.out && rm counts.out  && rm tmp.txt 

