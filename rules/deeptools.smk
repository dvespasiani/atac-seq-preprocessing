##=========================================
## 8. Run deeTools2 to get: 
## a) BAM files summary 
## b) correlation between samples
## c) fingerprint
## d) BAM coverage
##=========================================
rule_name = 'deeptools/'
dedup_bam = outdir + "post-alignment/{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bam"
index_bam = outdir + "post-alignment/{sample}-nochrM-encodefiltered-fixmate-rmorphanread-nodup.bai"

## First remove blacklisted regions and index resulting bam
rule deeptools_noblacklist:
  input:
    bam = dedup_bam,
    bai = index_bam
  output:
    outdir + rule_name + "{sample}-noblacklist.bam"
  params: 
   blacklist = blacklist
  group: 
   qc
  log:
    logs + rule_name + "{sample}-deeptools-noblacklist.log"
  shell:
   "bedtools intersect -nonamecheck -v -a {input.bam} -b {params.blacklist} > {output} 2> {log}"

rule deeptools_noblacklist_index:
  input:
   rules.deeptools_noblacklist.output
  output:
    outdir + rule_name + "{sample}-noblacklist.bai"
  group: 
   qc
  log:
    logs + rule_name + "{sample}-deeptools-noblacklist-index.log"
  shell:
   "samtools index -c {input} {output} 2> {log}"

rule deeptools_coverage:
  input:
    bam = dedup_bam,
    bai = index_bam
  output:
    outdir + rule_name + "{sample}-SeqDepthNorm.bw"
  group: 
   qc
  params:
    genome_size = genome_size
  log:
    logs + rule_name + "{sample}-deeptools-coverage.log"
  shell:
   """
   bamCoverage --bam {input.bam} --normalizeUsing RPGC \
   --effectiveGenomeSize {params.genome_size} --extendReads -o {output} 2> {log}
   """

rule deeptools_plot_coverage:
  input:
    bam = expand(dedup_bam,sample=sample),
    bai = expand(index_bam,sample=sample),
  output:
   outdir + rule_name + "samples-bam-coverage.png"
  params:
   read_minQ = read_minQ,
   sample = 25000000
  group: 
   qc
  log:
   logs + rule_name + "deeptools-plot-bam-coverage.log"
  shell:
   """
   plotCoverage --bamfiles {input.bam}\
   --smartLabels --numberOfSamples {params.sample} \
   --minMappingQuality {params.read_minQ} -o {output} 2> {log}
    """   

rule deeptools_fingerprint:
  input:
    bam = expand(dedup_bam,sample=sample),
    bai = expand(index_bam,sample=sample),
  output:
    fig = outdir + rule_name + "samples-plot-fingerprint.png",
    metrics = outdir + rule_name + "multiBAM-fingerprint-metrics.txt",
    rawcounts = outdir + rule_name + "multiBAM-fingerprint-rawcounts.txt"
  params:
   read_minQ = read_minQ
  group: 
   qc
  log:
   logs + rule_name + "deeptools-plot-fingerprint.log"
  shell:
    """
    plotFingerprint -b {input.bam} --plotFile {output.fig} \
    --outQualityMetrics {output.metrics} --outRawCounts {output.rawcounts} \
    --smartLabels --minMappingQuality {params.read_minQ} --numberOfProcessors 'max/2' 2> {log}
    """

rule computeGCbias:
  input:
    bam = dedup_bam,
    bai = index_bam
  output:
   content = outdir + rule_name + "{sample}-GC-content.txt",
   plot = outdir + rule_name + "{sample}-plot-GC-content.png"
  params:
   genome_size = genome_size,
   genome2bit = genome2bit_index,
   plot_format = 'png',
   threads = 8
  group: 
   qc
  log:
   logs + rule_name + "{sample}-deeptools-GC-content.log"
  shell:
   """
   computeGCBias -b {input.bam} --genome {params.genome2bit} \
   --effectiveGenomeSize {params.genome_size}  --plotFileFormat {params.plot_format} \
   --numberOfProcessors {params.threads} -o {output.content} --biasPlot {output.plot} 2> {log}
    """ 

# keep these inputs, as this rule will need all aligned/filtered bams
rule deeptools_summary:
  input:
    bam = expand(dedup_bam,sample=sample),
    bai = expand(index_bam,sample=sample),
  output:
    summary = outdir + rule_name + "multibam-summary.npz",
    readcounts = outdir + rule_name + "multibam-readcounts.txt"
  threads: 5
  group: 
   qc
  log:
    logs + rule_name + "deeptools-summary.log"
  shell:
   """
   multiBamSummary bins -p {threads} -b {input.bam} \
    -out {output.summary} --outRawCounts {output.readcounts} 2> {log}
    """

rule deeptools_correlation:
  input:
   rules.deeptools_summary.output.summary
  output:
    fig = outdir + rule_name  + "pearson-corr-multibam.png",
    matrix = outdir + rule_name  + "pearson-corr-multibamsum-matrix.txt"
  group: 
   qc
  log:
    logs + rule_name + "deeptools-correlation.log"
  shell:
    """
    plotCorrelation --corData {input} --plotFile {output.fig} \
    --outFileCorMatrix {output.matrix} --corMethod pearson \
    --whatToPlot heatmap --skipZeros --plotNumbers --colorMap RdYlBu 2> {log}
    """
