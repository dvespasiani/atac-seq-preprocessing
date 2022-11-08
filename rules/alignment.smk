## =================================================
##  1. Align to reference genome using Bowtie2
##  then sort bam file by genomic coordinates
## =================================================
priority = 1
rule align:
  input:
    r1 = rules.trim_adapter.output.r1,
    r2 = rules.trim_adapter.output.r2 
  output:
    outdir +  rulename_alignment + "{sample}.bam"
  params:
    index = genome_index
  log:
    logs + rulename_alignment + "{sample}.log"
  group:
    main
  priority:
    priority
  shell:
    """
    bowtie2 -q -X 2000 --mm -x {params.index} \
    -1 {input.r1} -2 {input.r2} 2> {log} \
    | samtools view -Su | samtools sort -n -o {output}
    """

## =====================================================================================##
##  Post-alignment QCsto remove unmapped, duplicated, multimapped and low quality reads ##
## =====================================================================================##

rule rmChrM:
  input:
    rules.align.output
  output:
    outdir +  rulename_alignment + "{sample}-nochrM.bam"
  group:
    main
  priority:
    priority
  log:
    logs + rulename_alignment + "{sample}-rmchrM.log"
  shell:
    "samtools view -h {input} | grep -v 'chrM' | samtools sort -o {output} 2>{log}"
 
rule encode_filters:
  input:
    rules.rmChrM.output
  output:
    outdir +  rulename_alignment + "{sample}-nochrM-encodefiltered.bam"
  group:
    main
  priority:
    priority
  params:
    read_minQ = read_minQ
  log:
    logs + rulename_alignment + "{sample}-encode-filters.log"
  shell:
    "samtools view -b -F 1804 -q {params.read_minQ} -f 2 -u {input} | samtools sort -n -o {output}  2> {log}"

rule fixmate:
  input:
    rules.encode_filters.output
  output:
    outdir +  rulename_alignment + "{sample}-nochrM-encodefiltered-fixmate.bam"
  group:
    main
  priority:
    priority
  log:
    logs + rulename_alignment  + "{sample}-fixmate.log"
  shell:
    "samtools fixmate -r {input} {output} 2> {log}"

rule rmOrphanread:
  input:
    rules.fixmate.output
  output:
    outdir + rulename_alignment + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread.bam"
  group:
    main
  priority:
    priority
  log:
    logs + rulename_alignment + "{sample}-rmOrphanread.log"
  shell:
    "samtools view -F 1804 -f 2 -u {input} | samtools sort -o {output} 2> {log}"

##============================
## Mark and remove duplicates
##============================
## remember to sort samfile by qname for this
## keep parameters USE_JDK_DEFLATER=true USE_JDK_INFLATER=true otherwise there could be java incompatibilities
rule markDups:
  input:
    rules.rmOrphanread.output 
  output:
    bam = outdir + rulename_alignment  + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread-dupmark.bam",
    dupQC = qcdir + rulename_alignment + "{sample}-duplicate-rate.qc"
  group:
    main
  priority:
    priority
  params:
   mem = "-Xmx4g"
  log:
    logs + rulename_alignment + "{sample}-markDups.log"
  shell:
    """
    MarkDuplicates I={input} O={output.bam} \
    METRICS_FILE={output.dupQC} VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> {log}
    """

rule dedup:
  input:
    rules.markDups.output.bam
  output:
    outdir + rulename_alignment + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bam"
  group:
    main
  priority:
    priority
  log:
    logs + rulename_alignment + "{sample}-dedup.log"
  shell:
    "samtools view -h -b  -F 1804 -f 2 {input} > {output} 2> {log}"

##============================
## 6.Index bam file
##============================
rule indexBam:
  input:
    rules.dedup.output
  output:
    outdir + rulename_alignment + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bai"
  group:
    main
  priority:
    priority
  log:
    logs + rulename_alignment + "{sample}-indexBam.log"
  shell:
    "samtools index -c {input} {output} 2> {log}"

## Tn5 shift reads
rule Tn5_shift:
  input:
    bam = rules.dedup.output,
    bai = rules.indexBam.output
  output:
    outdir + rulename_alignment + "{sample}-tn5-shifted.bam"
  group:
    main
  priority:
    priority
  log:
    logs + rulename_alignment + "{sample}-tn5-shifted.log"
  shell:
    "alignmentSieve -b {input.bam} -o {output} --ATACshift 2> {log}"

## pool together bams and create an extra big bam file and process this one along with the others.
## you'll need it for extracting the set of consensus peaks across all your samples
rule poolbams:
  input:
    samples = expand(outdir + rulename_alignment + "{sample}-tn5-shifted.bam",sample = sample)
  output:
    outdir + rulename_alignment + "{all_samples}-tn5-shifted.bam"
  log:
    logs + rulename_alignment + "{all_samples}-poolbams.log"
  group:
    main
  priority:
    priority
  shell:
    "samtools merge {output} {input.samples}"

rule Tn5_shifted_sort:
  input:
    outdir + rulename_alignment + "{combined_sample}-tn5-shifted.bam"
  output:
    outdir + rulename_alignment + "{combined_sample}-tn5-shifted-sorted.bam"
  group:
    main
  priority:
    priority
  log:
    logs + rulename_alignment + "{combined_sample}-tn5-shifted-sorted.log"
  shell:
    "samtools sort {input} -O BAM -o {output} 2> {log}"

rule index_Tn5Bams:
  input:
    rules.Tn5_shifted_sort.output
  output:
    outdir + rulename_alignment + "{combined_sample}-tn5-shifted-sorted.bam.bai"
  group:
    main
  priority:
    priority
  log:
    logs + rulename_alignment + "{combined_sample}-tn5-shifted-sorted-index.log"
  shell:
    "samtools index -c {input} {output} 2> {log}"