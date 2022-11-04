rule trimmomatic:
	input:
		r1 = fastqdir + "{sample}_R1_001.fastq.gz", 
		r2 = fastqdir + "{sample}_R2_001.fastq.gz", 
	output:
		r1 = outdir + rulename_trim + "{sample}-1-trimmed.fastq.gz",
		r2= outdir + rulename_trim + "{sample}-2-trimmed.fastq.gz",
		# reads where trimming entirely removed the mate
		r1_unpaired= outdir + rulename_trim + "{sample}-1-unpaired.fastq.gz",
		r2_unpaired= outdir + rulename_trim + "{sample}-2-unpaired.fastq.gz"
	params:
		adapters = adapters
	log:
		logs + rulename_trim + "{sample}.log"
	shell:
		"""
		trimmomatic  PE -phred33 \
		{input.r1} {input.r2} \
		{output.r1} {output.r1_unpaired} \
		{output.r2} {output.r2_unpaired} \
		ILLUMINACLIP:{params.adapters}:2:30:10 \
		LEADING:20 TRAILING:20 MINLEN:20 -trimlog {log}
		"""
