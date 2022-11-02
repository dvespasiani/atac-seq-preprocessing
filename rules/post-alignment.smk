## ===========================================================================================
##  Post-alignment QCsto remove unmapped, duplicated, multimapped and low quality reads
##  Removed also multimapping reads (MAPQ <=30) 
## ===========================================================================================
rule_name = 'post-alignment/'

rule rmChrM:
  input:
    rules.align.output
  output:
    outdir +  rule_name + "{sample}-nochrM.bam"
  group: 
    main
  log:
    logs + rule_name + "{sample}-rmchrM.log"
  shell:
    "samtools view -h {input} | grep -v 'chrM' | samtools sort -o {output} 2>{log}"
 
rule encode_filters:
  input:
    rules.rmChrM.output
  output:
    outdir +  rule_name + "{sample}-nochrM-encodefiltered.bam"
  group: 
    main
  params:
    read_minQ = read_minQ
  log:
    logs + rule_name + "{sample}-encode-filters.log"
  shell:
    "samtools view -b -F 1804 -q {params.read_minQ} -f 2 -u {input} | samtools sort -n -o {output}  2> {log}"

rule fixmate:
  input:
    rules.encode_filters.output
  output:
    outdir +  rule_name + "{sample}-nochrM-encodefiltered-fixmate.bam"
  group: 
    main
  log:
    logs + rule_name  + "{sample}-fixmate.log"
  shell:
    "samtools fixmate -r {input} {output} 2> {log}"

rule rmOrphanread:
  input:
    rules.fixmate.output
  output:
    outdir + rule_name + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread.bam"
  group: 
    main
  log:
    logs + rule_name  + "{sample}-rmOrphanread.log"
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
    bam = outdir + rule_name  + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread-dupmark.bam",
    dupQC = qcdir + "{sample}-duplicate-rate.qc"
  group: 
    main
  params:
   mem = "-Xmx4g"
  log:
    logs + rule_name + "{sample}-markDups.log"
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
    outdir + rule_name + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bam"
  group: 
    main
  log:
    logs + rule_name + "{sample}-dedup.log"
  shell:
    "samtools view -h -b  -F 1804 -f 2 {input} > {output} 2> {log}"

##============================
## 6.Index bam file
##============================
rule indexBam:
  input:
    rules.dedup.output
  output:
    outdir + rule_name + "{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bai"
  group: 
    main
  log:
    logs + rule_name + "{sample}-indexBam.log"
  shell:
    "samtools index -c {input} {output} 2> {log}"

## Tn5 shift reads
rule Tn5_shift:
  input:
    bam = rules.dedup.output,
    bai = rules.indexBam.output
  output:
    outdir + rule_name + "{sample}-tn5-shifted.bam"
  group: 
    main
  log:
    logs + rule_name + "{sample}-tn5-shifted.log"
  shell:
    "alignmentSieve -b {input.bam} -o {output} --ATACshift 2> {log}"

# rule poolbams:
#   input:
#     samples = expand(aligment_outdir + "{sample}_tn5_shifted.bam",sample = sample)
#   output:
#     aligment_outdir + "{merged_sample}_tn5_shifted.bam"
#   shell:
#     "samtools merge {output} {input.samples}"

rule Tn5_shifted_sort:
  input:
    rules.Tn5_shift.output
  output:
    outdir + rule_name + "{sample}-tn5-shifted-sorted.bam"
  group: 
    main
  log:
    logs + rule_name + "{sample}-tn5-shifted-sorted.log"
  shell:
    "samtools sort {input} -O BAM -o {output} 2> {log}"

rule index_Tn5Bams:
  input:
    rules.Tn5_shifted_sort.output
  output:
    outdir + rule_name + "{sample}-tn5-shifted-sorted.bam.bai"
  group: 
    main
  log:
    logs + rule_name + "{sample}-tn5-shifted-sorted-index.log"
  shell:
    "samtools index -c {input} {output} 2> {log}"