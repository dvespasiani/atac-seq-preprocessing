priority = 1
rule fastqc:
  input:
    inputdir + "{samples}_R1_001.fastq.gz", 
    inputdir + "{samples}_R2_001.fastq.gz"
  output:
    outdir + rulename_fastqc + "{samples}_R1_001_fastqc.zip",
    outdir +  rulename_fastqc + "{samples}_R2_001_fastqc.zip",
  log:
    logs +  rulename_fastqc + "{samples}.log"
  priority:
    priority
  group:
    main
  shell:
    "fastqc {input} --outdir=out/fastqc 2> {log}"
