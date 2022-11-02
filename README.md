# Table of Contents
- [Table of Contents](#table-of-contents)
  - [Project description](#project-description)
  - [Project set up](#project-set-up)
    - [Testing](#testing)
  - [Pipeline overview](#pipeline-overview)
    - [FastQC](#fastqc)
  - [Adaptor trimming](#adaptor-trimming)
  - [Genome alignment](#genome-alignment)
  - [Peak call](#peak-call)
  - [QCs](#qcs)
    - [Get peak and read counts](#get-peak-and-read-counts)
    - [Estimate library complexity](#estimate-library-complexity)
    - [Fraction Reads in Peaks (FRiP)](#fraction-reads-in-peaks-frip)
    - [TSS enrichment](#tss-enrichment)
    - [BAM Summary and Coverage (FROM HERE IT NEEDS TO BE UPDATED)](#bam-summary-and-coverage-from-here-it-needs-to-be-updated)
    - [GC bias](#gc-bias)
    - [Cumulative enrichment (BAM fingerprint)](#cumulative-enrichment-bam-fingerprint)
 
## Project description
This repo contains a snakemake pipeline to preprocess fastq files generated from bulk ATAC-seq experiments. Preprocessing steps mainly come from the [ENCODE ATAC-seq processing standards](https://www.encodeproject.org/atac-seq/). For full protocol specifications [check this google doc](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit). I have noticed they have changed it since last time (mainly polished it), so keep an eye on this. In addition to what reported by the ENCODE, I have also included some extra scripts to perfom and visualise quality control metrics.
## Project set up
To set up this pipeline you need to:
1. Modify the entries in the `config/snakemake-config.yaml` file, such as:
```
species: your-species
genome: your-species-genome
samples: 
    - your-sample1
    - your-sample2
basedir: your-basedir
etc....
```

**Importantly:** some softwares (e.g., MACS2, deeptools) require you to specify the effective genome size based on your sequencing read length for the species you are working with. You can either find this information [at this website](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html) or, alternatively, by running the `./bin/unique-kmers.py` script as:
```
python ./bin/unique-kmers.py -k <your-read-length> <path/to/genome/fasta/file.fa>
```
This script will return you the total estimated number of k-mers found in your species genome assembly. If you need further info on this script look at [MR Crusoe *et al.*, 2015](http://dx.doi.org/10.12688/f1000research.6924.1).

Be sure to have `pybedtools` python module installed in your python3. For the moment I have installed it in my own `PYTHONPATH` dir (which is also specified in my `~/.bash_profile`) and I have specified it in the `./bin/calculate-frip.py` script using the sys module in python. **This is not ideal**, but for the moment it works. I need to change it.

### Testing
If you want to test this pipeline on a subset of your fastq files to see if everything runs smoothly and all the desired results are produced then run `utils/subsample-files.sh -i <indir> -o <outdir> -n <numbreads>`. This script will extract the first `n = number of reads` from all the files listed in the `i = input dir`.

## Pipeline overview
The pipeline contained in this repo goes through the following steps:
1. FastQC
2. Adaptor trimming (Trimmomatic)
3. Alignment (Bowtie2)
4. Post-alignment filtering 
5. Peak-calling (MACS2)
7. QCs generation, e.g. using deepTools and other metrics

### FastQC
When running FastQC, you can expect 3 modules returning a warining/failure signal:
1. Per base sequence content because Tn5 has a strong sequence bias at the insertion site. 
2. Sequence Duplication Levels caused by PCR duplicates
3. Overrepresented sequences

## Adaptor trimming
Contrarily to what ENCODE did, for this step I used Trimmomatic to remove Illumina Nextera adapter sequeces. Fasta files containing the actual adapter sequences were obtained from  [here](https://github.com/timflutre/trimmomatic/blob/master/adapters/NexteraPE-PE.fa).<br/>
To identify adapters I've allowed for 2 max mismatches, a threshold of Q = 30 for PE palindrome read alignment (this can control for short adapter retentions at the 3' end of each read) and a threshold of Q = 10 for a simple alignment match between adapters and read. <br />
PS: In trimmomatic each match increase Q score by 0.6 whereas mismatches reduce Q by Q/10. Thus for a Q of 30 there must be around 50 matches and for a Q of 10 there are around 16 matches. <br />
Finally, I've removed all initial/terminal sequences having a phred score <20 and all reads that, after these quality steps were <20bp long. <br/>
Trimmomatic then splits fastq files in `*.paired.fastq` and `*.unpaired.fastq`. All downstream analyses are carried out only on paired trimmed reads.

## Genome alignment 
I used Bowtie2 to align reads to hg38 and panTro5 genomes. Differences with ENCODE alignment:
* I have specified a max fragment length of 5kb (encode sets it to 2kb)
* I have allowed dovetail (i.e. mates extending past each other)
* The maximum number of distinct, valid alignments was set to 1, so no multimapping was allowed (encode offers both solutions). Anyhow, Bowtie2 by default reports the best alignment, and in case of multimapping with same score a random selection is performed.
* I have ran Bowtie2 on a very sensitive mode (i.e., `--very-sensitive`) which by default means:
      * No mismatches are allowed 
      * The length of the seed substrings to align during multiseed alignment is 20 (ps: smaller values make alignment more sensitive)
      * Up to 20 consecutive seed extension attempts that can fail before Bowtie2 moves on, using the best alignment found (nb: a seed extension fails if it does not yield a new best or a new second-best alignment)
      * A maximum of 3 re-seeds of reads with repetitive seeds is allowed (i.e. Bowtie2 chooses a new set of reads with same length and same number of mismatches at different offsets and searches for more alignments).
See [Bowtie2 manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-i) for a more detailed description. <br/>
**NB**: Concordant alignment occurs when reads aligned with their expected mate orientation and with the expected range of distance between mates. Discordant alignment instead occurs when both mates have unique alignments, but the mates aren't in their expected relative orientations or aren't within their expected range ditance (or both).
## Peak call
Peak calling was performed using MACS2 on Tn5-shifted BAM files (from `ATACseqQC.R` script). Because coverage is very low for each sample I first merged all the technical replicates and then called the peaks on each sample (i.e. biological replicate).  <br/>
MACS2 flags set according to ENCODE guidelines.
Output files:
* $_default_peaks.narrowPeaks
BED6+4 file format. 
   * 5th column (i.e. score) is calculated as `int(-10*log10pvalue)`<br/>
   * 7th colum reports fold-change at peak summit <br/>
   * 8th and 9th colums report -log10pvalue and -log10qvalue respectively at peak summit <br/>
   * 10th column represents the relative summit position to peak start <br/>
* $_default_peaks.xlsx
Spreadsheet containing info about called peaks (i.e. chr; start; end; peak summit position; p/q-vals etc...)
* $_default_summits.bed
BED file containing the peak summits locations for every peak
* $_default_treat_pileup.bdg
bedGraph file containing the pileup signals of the treatment sample
* $_control_lambda.bdg
bedGraph file containing the local biases estimated for each genomic location from the control sample, or from treatment sample when the control sample is absent
<p>
These latter 2 files are then compared using `bdgcmp` subcommand to generate a signal track bedGraph file containing scores: p-value, q-value, log-likelihood, and log fold changes for each peak.

Chimp `chrom.sizes` file for bedtools was generated using:
```
samtools faidx panTro5.fa
cut -f 1,2 panTro5.fa.fai > chrom.sizes
```
whereas to get chrom.sizes for hg38 chromosomes simply run:
` wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes`


## QCs
Here I am listing a series of QCs the snakemake pipeline will automatically run for you plus some others you can perform by running some executables scripts within the `bin/` directory. Some of these QCs come from the deeptools suite of commands, check out the [deepTools manual](https://deeptools.readthedocs.io/en/develop/index.html) for a complete understanding. <br/>
### Get peak and read counts
Use `./bin/get-counts.sh -i <input_dir> -o <out_file>`  to obtain number of reads/peaks for each `*.bam$` and `*narrowPeak` file within the respective `out/preprocessing/` directories.
### Estimate library complexity
To get metrics such as PBC1/PBC2 and NRF (see ) run:
```
./bin/estimate-lib-complexity.sh -i <input_dir> -o <out_dir>
```
This will create a tab-separated `library-complexity.txt` file in your specified output dir with all those information.
### Fraction Reads in Peaks (FRiP)
It represents the proportion of all mapped reads that fall into the called peak regions. FRiP scores positively correlates with the number of regions. According to ENCODE: "*FRiP should be >0.3, though values greater than 0.2 are acceptable*". <br/>
FRiP calculation is already included in the snakemake pipeline, see `rule_frip` in `rules/peak-calling.smk`. However, to generate a barchart with the FRiP scores for each of your samples run:
```
bin/plot-frip-summary.R -i input/dir/containing/frip/result/txt/files
```
### TSS enrichment
This metric is used as another signal to noise indication. Reads around a reference set of TSSs are collected to form an aggregate distribution of reads centered on the TSSs and extending to 1000 bp in either direction (for a total of 2000bp). This distribution is then normalized by taking the average read depth in the 100 bps at each of the end flanks of the distribution (for a total of 200bp of averaged data) and calculating a fold change at each position over that average read depth. This means that the flanks should start at 1, and if there is high read signal at transcription start sites (highly open regions of the genome)there should be an increase in signal up to a peak in the middle. We take the signal value at the center of the distribution after this normalization as our TSS enrichment metric.

### BAM Summary and Coverage (FROM HERE IT NEEDS TO BE UPDATED)
Following BAM filtering I removed reads overlapping ENCODE blacklisted regions and then calculated the bam coverage with `bamCoverage` for all samples. Coverage is calculated as the number of reads per bin (i.e. short consecutive counting windows of defined size) and it can be scaled in RPKM,CPM. <br/>
Resulting summary was plotted using `plotCoverage` to visually assess the sequencing depth of each sample after sampling 25*10^6 reads. This command counts the number of overlapping reads and returns 2 plots indicating:
1. The frequencies of the observed read coverages per sample
2. The fraction of the genome covered by >= a given number of reads 

To investigate the correlation between the samples read coverages I've used `multiBamSummary`, which computes the read coverages over the entire genome and/or specific genomic regions for >=2 BAM files. Resulting pairwise Pearson correlation values are plotted into a heatmap  with `plotCorrelation`. 

### GC bias
This QC step tests the assumption of an expected uniform distribution of sequenced reads across the genome, regardless of their base-pair composition. However, PCR steps performed during library preparation can enrich for GC-rich fragments and this will influence the sequencing results.<br/>
The deepTools command `computeGCbias` first calculates the expected GC profile from the given genome by calculating the GC % for DNA fragments of a fixed size. It then compares the expectations with the observed GC content within the sequenced reads. <br/> 
The resulting plots report the number of reads overlapping genomic regions of increasing GC content and the log2(observed/expected) ratio of reads per GC content region. If there is no GC bias then the observed and expected GC profiles will be similar, i.e. ratio of observed/expected GC content per fragment length would = 1. <br/>
In this experiment there is a linear increase in the exp/obs read coverage across all samples. Possible solutions:
1. Run a simple, non-conservative peak calling on the uncorrected BAM file first to obtain a BED file of peak regions that are then supplied to `computeGCbias` to limit the calculation
2. Run the `correctGCbias` deepTool, which removes reads from regions of too high-coverage and add reads to regions of low-coverage
2. Leave it as it is since the bias is consistent across all samples and it might be a biological effect 

### Cumulative enrichment (BAM fingerprint)
<p> The deepTools `plotFingerprint` samples indexed BAM files and calculates the samples cumulative enrichment by counting all reads overlapping a window (bin) of a given length. The resulting plot is useful to assess how well the signal of the sequenced reads can be differentiated from the background distribution of reads in the control input sample. <br/>

An ideal input with perfect uniform distribution of reads along the genome (i.e. without enrichments in specific open chromatin regions) and infinite sequencing coverage would generate a straight diagonal line. On the other hand, a strong enrichment in very specific (i.e. fragment size) chromatin regions will result in a prominent and steep rise of the cumulative sum towards the highest rank. This means that a big chunk of reads from the test (ChIP/ATAC) sample is located in few bins which corresponds to high, narrow enrichments typically seen for transcription factors. <br/>

**PS**: check where cumulative curve starts on the plots as this will give you an indication of the percentage of the genome that was not sequenced at all (i.e. bins containing 0 reads). However, such percentage can be quite high, especially for extremely high enrichment which will result in the vast majority of reads occurring within few peaks, [see examples](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html).  </p>
