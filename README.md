# Table of Contents
- [Table of Contents](#table-of-contents)
  - [Project description](#project-description)
  - [Project set up](#project-set-up)
  - [Running the pipeline](#running-the-pipeline)
    - [Interactively](#interactively)
    - [In the background](#in-the-background)
  - [Pipeline overview](#pipeline-overview)
    - [FastQC](#fastqc)
    - [Adaptor trimming](#adaptor-trimming)
    - [Genome alignment](#genome-alignment)
    - [Peak calling](#peak-calling)
  - [Quality controls](#quality-controls)
    - [Get peak and read counts](#get-peak-and-read-counts)
    - [Estimate library complexity](#estimate-library-complexity)
    - [Get summary alignment results](#get-summary-alignment-results)
    - [Peak general QCs](#peak-general-qcs)
    - [Fraction Reads in Peaks (FRiP)](#fraction-reads-in-peaks-frip)
    - [TSS enrichment](#tss-enrichment)
    - [BAM Summary and Coverage](#bam-summary-and-coverage)
    - [GC bias](#gc-bias)
    - [Cumulative enrichment (BAM fingerprint)](#cumulative-enrichment-bam-fingerprint)
 

----
## Project description
This repo contains a snakemake pipeline to preprocess fastq files generated from bulk ATAC-seq experiments. Preprocessing steps mainly come from the [ENCODE ATAC-seq processing standards](https://www.encodeproject.org/atac-seq/). For full protocol specifications [check this google doc](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit). I have noticed they have changed it since last time (mainly polished it), so keep an eye on this. In addition to what reported by the ENCODE, I have also included some extra scripts to perfom and visualise quality control metrics.
## Project set up
To set up this pipeline you need to:
1. Modify the entries in the `config/snakemake-config.yaml` file to your needs, such as:
```
species: your-species
genome: your-species-genome
samples: 
    - your-sample1
    - your-sample2
basedir: your-basedir
fastqdir: your-fastqdir-containing-fastq-files
etc....

```
2. Make sure you have all the information for your species, such as:
   * a chrom.sizes file containing the chromosome sizes. This can be obtained either from UCSC (the link should be something like `http://hgdownload.soe.ucsc.edu/goldenPath/<your-species-assembly>/bigZips/<your-species-assembly>.chrom.sizes`) or, alternatively, by running
```
 samtools faidx <your-species-assembly>.fa
 cut -f 1,2 <your-species-assembly>.fa.fai > <your-species-assembly>.chrom.sizes
 ``` 
   * the effective genome size for your species of interest which is based on the length of your sequencing reads. Again, you can find this info either [at this website](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html) or, alternatively, by running
```
 python ./bin/unique-kmers.py -k <your-read-length> <path/to/genome/fasta/file.fa>
 ``` 
This script will return you the total estimated number of k-mers found in your species genome assembly. If you need further info on this script look at [MR Crusoe *et al.*, 2015](http://dx.doi.org/10.12688/f1000research.6924.1). <br/>

3. Make sure you have all the softwares installed in your R/python envs:
   * The `pybedtools` python module. For the moment I have installed it in my own `PYTHONPATH` dir (which is also specified in my `~/.bash_profile`) and I have specified this path in the `./bin/calculate-frip.py` script using the sys module in python. **However, this is not ideal**, but for the moment it works. I need to change it.
   * R libraries such as `yaml`, `data.table`, `argparse`, `GenomicAlignments`, `ATACseqQC`, `csaw`, `GenomicFeatures`, `TxDb.<your-species>.UCSC.<your-species-assembly>.knownGene`. All the other libraries should be pretty standards

4. Once you know how many samples you are preprocessing, replace the following lines within the `utils/r-utils.R`:
```
qualitative_palette = brewer.pal.info[brewer.pal.info$category == 'qual',]
sample_palette = sample(unlist(mapply(brewer.pal, qualitative_palette$maxcolors, rownames(qualitative_palette))),length(samples))
names(sample_palette) = samples
```
with something like:
```
sample_palette =  c('your','palette')
```

## Running the pipeline

### Interactively
To test this pipeline on a subset of your fastq files to see if everything runs smoothly and all the desired results are produced, you need to run:

```
utils/subsample-files.sh -i <indir> -o <outdir> -n <numbreads>
```

This script will extract the first `n = number of reads` from all the files listed in the `i = input dir`. You can then run the entire pipeline interactively on this subset of reads. However, if you do so, remember to correctly specify the location of your subsampled fastq files in the `config/snakemake-config.yaml` file (hint: it's the `fastqdir` directive).
Next, move the `runInteractively.sh` script outside this project directory and set up your interactive session as:

```
chmod +x runInteractively.sh
./runInteractively.sh -p <your-project-name> # to change slurm salloc arguments (e.g., mem/time etc..) go to config/cluster_config.yaml
eval "$(cat "$modules")"  # this loads all the modules listed in the config/modules.txt file
```
Then simply run:
```
snakemake --cores 8 -s snakefile-preprocess.smk
```
and you should get your results. <br/>

### In the background

----
## Pipeline overview
The pipeline contained in this repo goes through the following steps:
1. FastQC
2. Adaptor trimming (Trimmomatic)
3. Alignment (Bowtie2)
4. Post-alignment filtering 
5. Peak-calling (MACS2)
7. QCs generation, e.g. using deepTools and other metrics
### FastQC
For ATAC-seq experiments, when running FastQC you can expect 3 modules returning a warining/failure signal:
1. Per base sequence content because Tn5 has a strong sequence bias at the insertion site. 
2. Sequence Duplication Levels caused by PCR duplicates
3. Overrepresented sequences
Check them anyhow in case there is something extremely weird in your dataset.
### Adaptor trimming
I am using Trimmomatic to remove Illumina Nextera adapter sequeces which I have obtained [here](https://github.com/timflutre/trimmomatic/blob/master/adapters/NexteraPE-PE.fa). If you have other adapter sequences then simply make a fasta file with those sequences and change the `adapters:` directive in the `config/snakemake-config.yaml` file. To identify and then remove adapters, I am allowing for 2 max mismatches, a threshold of Q = 30 for PE palindrome read alignment (this can control for short adapter retentions at the 3' end of each read) and a threshold of Q = 10 for a simple alignment match between adapters and read. Finally, I am removing all initial/terminal sequences having a phred score <20 and all reads that, after these quality steps are < 20 bp long. <br/>
**PS:** from the trimmomatic manual, each match increases the Q score by 0.6 whereas mismatches reduces it by Q/10. Thus when Q = 30, there should be around 50 matches between your sequence and the adapters, whereas Q = 10, there should be around 16 matches. 
### Genome alignment 
I am using Bowtie2 to align reads to the genome assembly of the species of interest. Bowtie2 alignment parameters are defined as per ENCODE.
### Peak calling
To call peaks of chromatin accessiblity I am using MACS2 on Tn5-shifted BAM files. MACS2 parameters are defined as per ENCODE.
## Quality controls
Here I am listing a series of QCs the snakemake pipeline will automatically run for you plus some others you can perform by running some executables scripts within the `bin/` directory. Some of these QCs come from the deeptools suite of commands, check out the [deepTools manual](https://deeptools.readthedocs.io/en/develop/index.html) for a complete understanding. <br/>
All the QCs I am generating extra are by default located in the `out/preprocessing/plot/atac-seq-qc` directory. You can change this if you want to.
### Get peak and read counts
Use `./bin/get-counts.sh -i <input_dir> -o <out_file>`  to obtain number of reads/peaks for each `*.bam$` and `*narrowPeak` file within the respective `out/preprocessing/` directories.
### Estimate library complexity
To get metrics such as PBC1/PBC2 and NRF ([see here](https://www.encodeproject.org/data-standards/terms/#library)) run:
```
./bin/estimate-lib-complexity.sh -i <input_dir> -o <out_dir>
```
This will create a tab-separated `library-complexity.txt` file in your specified output dir with all those information.
### Get summary alignment results
To plot the Bowtie2 alignment results (i.e., the % of reads aligned) for each sample, run:
```
Rscript ./bin/plot-alignment-summary.R -i path/to/logs/alignment/dir/
```
**NB:** The files you should access to, by default, are located in the `logs/alignment` directory as per snakemake.

### Peak general QCs
To plot some other general QCs for your set of peaks run:
```
Rscript ./bin/plot-peak-qcs.R 
```
At the moment this script will only return you a plot with the distribution of peak sizes for your set of peaks. If I will think/come across with other QCs I will incorporate them in here
### Fraction Reads in Peaks (FRiP)
It represents the proportion of all mapped reads that fall into the called peak regions. FRiP scores positively correlates with the number of regions. According to ENCODE: "*FRiP should be >0.3, though values greater than 0.2 are acceptable*". <br/>
FRiP calculation is already included in the snakemake pipeline, see `rule_frip` in `rules/peak-calling.smk`. However, to generate a barchart with the FRiP scores for each of your samples run:

```
Rscript ./bin/plot-frip-summary.R -i path/to/input/dir/with/frip/result/
```
### TSS enrichment
To calculate and plot the TSS enrichment run:
```
Rscript ./bin/plot-tss-enrich.R
```
### BAM Summary and Coverage
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
