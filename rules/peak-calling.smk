##===========================##
## 10. Peak calling (MACS2)  ##
##===========================##

rule peakCalling:
  input:
    bam = rules.Tn5_shifted_sort.output,
    bai = rules.index_Tn5Bams.output
  output:
    peaks_xls = outdir + rulename_peak + "{combined_sample}-macs2_peaks.xls",
    pileup = outdir + rulename_peak + "{combined_sample}-macs2_treat_pileup.bdg",
    lamb = outdir + rulename_peak + "{combined_sample}-macs2_control_lambda.bdg",
    summit = outdir + rulename_peak + "{combined_sample}-macs2_summits.bed",
    narrowPeak = outdir + rulename_peak + "{combined_sample}-macs2_peaks.narrowPeak"
  group: 
    main
  params:
    name = "{combined_sample}-macs2",
    outdir = outdir + rulename_peak,
    fragment_size = fragment_size,
    shift = shift,
    genome_size = genome_size,
    pval_thresh = pval_thresh
  log:
    logs + rulename_peak + "{combined_sample}-peak-calling.log"
  shell:
    """
    macs2 callpeak --format BAMPE --treatment {input.bam} \
    --keep-dup all --outdir {params.outdir} --name {params.name} --shift {params.shift} \
    --nomodel -B --SPMR --extsize {params.fragment_size} \
    --pvalue {params.pval_thresh} --call-summits -g {params.genome_size} 2> {log}
    """

rule rm_blacklisted_peaks:
  input:
    rules.peakCalling.output.narrowPeak
  output:
    outdir + rulename_peak + "{combined_sample}-macs2-peaks-filtered.narrowPeak.gz"
  params:
    blacklist = blacklist
  group: 
    main
  log:
    logs + rulename_peak + "{combined_sample}_blacklist_removed.log"
  shell:
    """
    bedtools intersect -v -a {input} -b {params.blacklist} \
    | awk 'BEGIN{{OFS="\\t"}} {{if ($5>1000) $5=1000; print $0}}' \
    | grep -P 'chr[\dXY]+[ \t]'  | gzip -nc > {output} 2> {log}
    """

## Sort NarrowPeaks by column 8 in descending order,
## replace long peak names in column 4 with peak_<peakRank> and keep only first NPEAKS
rule sort_peaks:
  input:
    rules.rm_blacklisted_peaks.output
  output:
    outdir + rulename_peak + "{combined_sample}-macs2-peaks-filtered-sorted.narrowPeak.gz"
  log:
     logs + rulename_peak + "{combined_sample}-sorted-narrowpeaks.log"
  params:
    npeaks=npeaks
  group: 
    main
  shell:
    """
    gunzip -nc {input} | sort -k 8gr,8gr | \
    awk 'BEGIN{{OFS="\\t"}}{{$4="peak_"NR ; print $0}}' \
    | head -n {params.npeaks} | gzip -nc > {output} 2> {log}
    """
  
##------------------------------------------------
## Get peak Fold Enrichments (FE) and 
## then convert it to bigWig format for display
##------------------------------------------------
rule FE_peak_signal_tracks:
  input:
   treatment = rules.peakCalling.output.pileup,
   control = rules.peakCalling.output.lamb
  output:
    outdir + rulename_peak + "{combined_sample}-macs2_FE.bdg"
  group: 
    main
  log:
   logs + rulename_peak + "{combined_sample}-fe-peak-signal-track.log"
  shell:
    """
    macs2 bdgcmp -t {input.treatment} -c {input.control} \
    --ofile {output} -m FE  2> {log}
    """

rule clean_FE_signal:
  input:
   rules.FE_peak_signal_tracks.output
  output:
    outdir + rulename_peak + "{combined_sample}-fe-signal.bedgraph"
  log:
   logs + rulename_peak + "{combined_sample}-clean-FE-peaks.log"
  params:
    chrom_sizes = chrom_sizes
  group: 
    main
  shell:
   """
   slopBed -i {input} -g {params.chrom_sizes} -b 0 \
   | bedClip stdin {params.chrom_sizes} {output}  2> {log}
   """

rule sort_FE_bedGraph:
  input:
   rules.clean_FE_signal.output 
  output:
    outdir + rulename_peak + "{combined_sample}-fe-signal-sorted.bedgraph"
  group: 
    main
  log:
   logs + rulename_peak + "{combined_sample}-sort-bedGraph.log"
  shell:
   "sort -k1,1 -k2,2n {input} > {output} 2> {log}"

rule FE_bedGraph2bigWig:
  input:
   rules.sort_FE_bedGraph.output
  output:
   outdir + rulename_peak + "{combined_sample}-fe-signal.bigwig"
  group: 
    main
  log:
   logs + rulename_peak + "{combined_sample}-bedGraph2bigWig.log"
  params:
    chrom_sizes = chrom_sizes
  shell:
   """
   bedGraphToBigWig {input} {params.chrom_sizes} {output} 2> {log}
   """

##------------------------------------
## now repeat these latter steps 
## but to calculate peak p-values
##-----------------------------------

## sval counts the number of tags per million in the (compressed) BED file
rule sval:
  input:
   rules.peakCalling.output.narrowPeak
  output:
   outdir + rulename_peak + "{combined_sample}-ppois-sval"
  group: 
    main
  log:
   logs + rulename_peak + "{combined_sample}-sval-calculation.log"
  shell:
   """
   (wc -l <(gunzip -nc {input}) \
   | awk '{{printf "%f", $1/1000000}}' >{output} ) 2> {log}
   """

rule ppois_peak_signal_tracks:
  input:
   sval = rules.sval.output,
   treatment = rules.peakCalling.output.pileup,
   control = rules.peakCalling.output.lamb
  output:
    outdir + rulename_peak + "{combined_sample}-macs2_ppois.bdg"
  group: 
    main
  log:
    logs + rulename_peak + "{combined_sample}-ppois-peak-signaltrack.log"
  run:
   for filename in glob.glob(outdir + rulename_peak + '*-ppois-sval'):
      sval = open(filename).read()
      shell("""macs2 bdgcmp -t {input.treatment} -c {input.control} --ofile {output} \
      -m ppois -S {sval} 2> {log}""")

rule clean_ppois_signal:
  input:
   rules.ppois_peak_signal_tracks.output
  output:
    outdir + rulename_peak + "{combined_sample}-ppois-signal.bedgraph"
  group: 
    main
  log:
   logs + rulename_peak + "{combined_sample}-clean-ppois-peaks.log"
  params:
    chrom_sizes = chrom_sizes
  shell:
   """
   slopBed -i {input} -g {params.chrom_sizes} -b 0 \
   | bedClip stdin {params.chrom_sizes} {output} 2> {log}
   """

rule sort_ppois_bedGraph:
  input:
   rules.clean_ppois_signal.output
  output:
    outdir + rulename_peak + "{combined_sample}-ppois-signal-sorted.bedgraph"
  group: 
    main
  log:
   logs + rulename_peak + "{combined_sample}-sort-ppois-bedGraph.log"
  shell:
   "sort -k1,1 -k2,2n {input} > {output} 2> {log}"

rule ppois_bedGraph2bigWig:
  input:
   rules.sort_ppois_bedGraph.output
  output:
   outdir + rulename_peak + "{combined_sample}-ppois-signal.bigwig"
  group: 
    main
  log:
   logs + rulename_peak + "{combined_sample}-ppois-bedGraph2bigWig.log"
  params:
    chrom_sizes = chrom_sizes 
  shell:
   """
   bedGraphToBigWig {input} {params.chrom_sizes} {output} 2> {log}
   """

##----------------------------------------------------
## FRiP on pooled bam files and filtered narrowPeaks
##----------------------------------------------------
rule frip:
  input:
   bams = rules.Tn5_shifted_sort.output,
   peaks = rules.sort_peaks.output
  output:
    qcdir + "{combined_sample}-frip.txt"
  group: 
    main
  log:
   logs + rulename_peak + "{combined_sample}-FRiP.log"
  run:
   shell("python3 ./bin/calculate-frip.py {input.bams} {input.peaks} > {output} 2> {log}")

# ##------------------------
# ## Call consensus peak
# ##------------------------
# rule consensus_peak:
#   input:
#    peak_merged_sample = peak_dir + "merged_macs2_default_peaks_filtered_sorted.narrowPeak.gz",
#    peak_other_samples = expand(peak_dir + "{augmented_samples}_macs2_default_peaks_filtered.narrowPeak.gz", sample=config["samples"])
#   output:
#    "output/PeakCalling/ConsensusPeaks/consensus_peak.bed.gz"
#   params:
#    min_overlap = 0.5,
#    sample_names = config["samples"]
#   group:
#     main
#   log:
#    "logs/PeakCalling/Consensus_peak.log"
#   shell:
#    """
#    bedtools intersect -f {params.min_overlap} -r \
#    -a {input.peak_merged_sample} -b {input.peak_other_samples} \
#    | gzip -nc > {output} 2> {log}
#    """
