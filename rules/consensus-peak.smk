priority = 1
rule consensus_peak:
    input:
        sample_peaks = expand(outdir + rulename_peak + "{sample}-macs2-peaks-filtered-sorted.narrowPeak.gz",sample = sample),
        combined_peaks = expand(outdir + rulename_peak + "{all_samples}-macs2-peaks-filtered-sorted.narrowPeak.gz",all_samples = all_samples)
    output:
        outplot = plots + rulename_qc + 'support-consensus-peak.pdf',
        outfile = outdir + rulename_conspeak + 'consensus-peak.bed.gz'
    log:
        logs + rulename_conspeak + "consensus-peak.log"
    group:
        main
    priority:
        priority
    script:
        "../bin/get-consensus-peaks.R"
