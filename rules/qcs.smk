rule count_bam_reads:
    input:
        outdir + rulename_postalign
    output:
        tables + 'number-reads-bam-files.txt'
    log:
        logs + rulename_postalign + "numb-reads-bams.log"
    group:
        qc
    shell:
        "./bin/get-counts.sh -i {input} -o {output}"


rule count_peaks:
    input:
        outdir + rulename_peak
    output:
        tables + 'number-peaks.txt'
    log:
        logs + rulename_peak + "number-peaks.log"
    group:
        qc
    shell:
        "./bin/get-counts.sh -i {input} -o {output}"


rule alignment_summary:
    input:
        logs + rulename_align
    output:
        plots + 'bowtie2-alignment-qc.pdf'
    log:
        logs + rulename_align + "alignment-summary.log"
    group:
        qc
    shell:
        "Rscript ./bin/plot-alignment-summary.R {input} {output}"


rule estim_lib_complex:
    input:
        outdir + rulename_postalign
    output:
        tables + 'library-complexity.txt'
    log:
        logs + rulename_postalign + "library-complexity.log"
    group:
        qc
    shell:
        "./bin/estimate-lib-complexity.sh -i {input} -o {output}"

rule frip_summary:
    input:
        qcdir
    output:
        plots + 'frip-summary-qc.pdf'
    log:
        logs + rulename_peak + "frip-summary-qc.log"
    group:
        qc
    shell:
        "Rscript ./bin/plot-frip-summary.R {input} {output}"
        

rule peak_qc:
    input:
        outdir + rulename_peak
    output:
        plots + 'extra-peak-qcs.pdf'
    log:
        logs + rulename_peak + "extra-peak-qcs.log"
    group:
        qc
    shell:
        "Rscript ./bin/plot-peak-qcs.R {input} {output}"
        

rule tss_enrich:
    input:
        outdir + rulename_postalign
    output:
        plots + 'tss-enrichment.pdf'
    log:
        logs + rulename_peak + "tss-enrichment.log"
    group:
        qc
    shell:
        "Rscript ./bin/plot-tss-enrich.R {input} {output}"
        

