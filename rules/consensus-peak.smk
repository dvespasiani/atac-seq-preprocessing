rule consensus_peak:
    input:
        outdir + rulename_peak
    output:
        outplot = plots + 'support-consensus-peak.pdf',
        outfile = outdir + rulename_conspeak + 'consensus-peak.bed.gz'
    log:
        logs + rulename_conspeak + "consensus-peak.log"
    group:
        qc
    shell:
        "Rscript ./bin/get-consensus-peaks.R {input} {output.outplot} {output.outfile}"
