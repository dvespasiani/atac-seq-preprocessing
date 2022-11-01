rule_name = 'fastqc/'

rule fastqc:
  input:
    fastqdir + "{samples}_R1_001.fastq.gz", 
    fastqdir + "{samples}_R2_001.fastq.gz"
  output:
    outdir + rule_name + "{samples}_R1_001_fastqc.zip",
    outdir + rule_name + "{samples}_R2_001_fastqc.zip",
  log:
    logs + rule_name + "{samples}.log"
  shell:
    "fastqc {input} --outdir=out/preprocessing/fastqc/ 2> {log}"
