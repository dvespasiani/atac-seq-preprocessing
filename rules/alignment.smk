## =================================================
##  1. Align to reference genome using Bowtie2
##  then sort bam file by genomic coordinates
## =================================================
rule align:
  input:
    r1 = rules.trimmomatic.output.r1,
    r2 = rules.trimmomatic.output.r2 
  output:
    outdir +  rulename_align + "{sample}.bam"
  params:
    index = genome_index
  group:
    main
  log:
    logs + rulename_align + "{sample}.log"
  shell:
    "bowtie2 -q -X 2000 --mm -x {params.index} \
    -1 {input.r1} -2 {input.r2} 2> {log} \
    | samtools view -Su | samtools sort -n -o {output}"

