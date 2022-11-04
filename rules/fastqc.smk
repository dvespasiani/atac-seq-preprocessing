rule fastqc:
  input:
    fastqdir + "{samples}_R1_001.fastq.gz", 
    fastqdir + "{samples}_R2_001.fastq.gz"
  output:
    outdir + rulename_fastqc + "{samples}_R1_001_fastqc.zip",
    outdir + rulename_fastqc + "{samples}_R2_001_fastqc.zip",
  log:
    logs + rulename_fastqc + "{samples}.log"
  shell:
    "fastqc {input} --outdir=out/preprocessing/fastqc/ 2> {log}"
