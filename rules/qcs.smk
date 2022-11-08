priority = 0 
rule count_bam_reads:
    input:
        outdir + rulename_alignment + "all_samples-tn5-shifted-sorted.bam",
    output:
        tables + rulename_qc + 'number-reads-bam-files.txt'
    log:
        logs + rulename_alignment + "numb-reads-bams.log"
    group:
        qc
    priority:
        priority
    script:
        "../bin/get-counts.py"
        
rule count_peaks:
    input:
        outdir + rulename_peak + "all_samples-macs2-peaks-filtered-sorted.narrowPeak.gz"
    output:
        tables + rulename_qc + 'number-peaks.txt'
    log:
        logs + rulename_peak + "number-peaks.log"
    group:
        qc
    priority:
        priority
    script:
        "../bin/get-counts.py"
        
rule alignment_summary:
    input:
        logs + rulename_alignment + "all_samples-tn5-shifted-sorted-index.log"
    output:
        plots + rulename_qc + "bowtie2-alignment-summary.pdf"
    log:
        logs + rulename_alignment + "alignment-summary.log"
    group:
        qc
    priority:
        priority
    script:
        "../bin/plot-alignment-summary.R"

rule estim_lib_complex:
    input:
        outdir + rulename_alignment + "all_samples-tn5-shifted-sorted.bam"
    output:
        tables + rulename_qc + 'library-complexity.txt'
    log:
        logs + rulename_alignment + "library-complexity.log"
    group:
        qc
    priority:
        priority
    script:
        "../bin/estimate-lib-complexity.py"

rule frip_summary:
    input:
        qcdir + "all_samples-frip.txt"
    output:
        plots + rulename_qc + 'frip-summary-qc.pdf'
    log:
        logs + rulename_peak + "frip-summary-qc.log"
    group:
        qc
    priority:
        priority
    script:
        "../bin/plot-frip-summary.R"
        
rule peak_qc:
    input:
        outdir + rulename_peak + "all_samples-macs2-peaks-filtered-sorted.narrowPeak.gz"
    output:
        plots + rulename_qc + 'extra-peak-qcs.pdf'
    log:
        logs + rulename_peak + "extra-peak-qcs.log"
    group:
        qc
    priority:
        priority
    script:
        "../bin/plot-peak-qcs.R"
        
rule tss_enrich:
    input:
        outdir + rulename_alignment + "all_samples-tn5-shifted-sorted.bam"
    output:
        plots + rulename_qc + 'tss-enrichment.pdf'
    log:
        logs + rulename_peak + "tss-enrichment.log"
    group:
        qc
    priority:
        priority
    script:
        "../bin/plot-tss-enrich.R"
        

