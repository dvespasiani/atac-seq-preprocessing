# ATAC-seq preprocessing pipeline
----
# Table of Contents
- [ATAC-seq preprocessing pipeline](#atac-seq-preprocessing-pipeline)
- [Table of Contents](#table-of-contents)
  - [Project description](#project-description)
  - [Project set up](#project-set-up)
  - [Running the pipeline](#running-the-pipeline)
    - [Interactively](#interactively)
    - [In the background](#in-the-background)
  - [**PS:**  Be sure to correctly specify the name AND location of your python virtual environment in the corresponding line in the `utils/run-on-cluster.sh` file! There should be no need to change anything else in there.](#ps--be-sure-to-correctly-specify-the-name-and-location-of-your-python-virtual-environment-in-the-corresponding-line-in-the-utilsrun-on-clustersh-file-there-should-be-no-need-to-change-anything-else-in-there)
  - [Pipeline overview](#pipeline-overview)
    - [FastQC](#fastqc)
    - [Adaptor trimming](#adaptor-trimming)
    - [Genome alignment](#genome-alignment)
    - [Peak calling](#peak-calling)
  - [Quality controls](#quality-controls)
      - [BAM Summary and Coverage](#bam-summary-and-coverage)
      - [GC bias](#gc-bias)
      - [Cumulative enrichment (BAM fingerprint)](#cumulative-enrichment-bam-fingerprint)
 

----
## Project description
This repo contains a snakemake pipeline to preprocess fastq files generated from bulk ATAC-seq experiments. Preprocessing steps mainly come from the [ENCODE ATAC-seq processing standards](https://www.encodeproject.org/atac-seq/). For full protocol specifications [check this google doc](https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit). I have noticed they have changed it since last time (mainly polished it), so keep an eye on this. In addition to what reported by the ENCODE, I have also included some extra scripts to perfom and visualise quality control metrics.

## Project set up
To set up this pipeline you need to:
* **FIRST AND FOREMOST** = have your tmp directory within the vast/scratch filesystem (if you dont have permissions there, contact IT)
If you do already have one then be sure the path it's in your `~/.bash_profile` and if not then save this line in there:
```
export TMPDIR=/vast/scratch/users/<your-username>/tmp
```
This is fundamental to have, if not for managing files and resourses in the cluster, simply because some commands of this pipeline will fail as they will fill up your base disk quota in no time. Anyhow, once you have it, modify the corresponding entry within the `snakemake-config.yaml` file (see right below here).

* Change the entries in the `config/snakemake-config.yaml` file to suit your needs, such as:
```
species: your-species
genome: your-species-genome
samples: 
    - your-sample1
    - your-sample2
basedir: your-basedir
fastqdir: your-dir-containing-fastq-files
etc....

```
Note that I am using full paths for most directives in the config yaml file. You could make them relative to a base directory to avoid repetition but you'll need to combine them in the main `snakemake-preprocess.smk` file.

* Make sure you have all the information for your species, such as:
   - a chrom.sizes file containing the chromosome sizes. This can be obtained either from UCSC (the link should be something like `http://hgdownload.soe.ucsc.edu/goldenPath/<your-species-assembly>/bigZips/<your-species-assembly>.chrom.sizes`) or, alternatively, by running
```
 samtools faidx <your-species-assembly>.fa
 cut -f 1,2 <your-species-assembly>.fa.fai > <your-species-assembly>.chrom.sizes
 ``` 
   - the effective genome size for your species of interest which is based on the length of your sequencing reads. Again, you can find this either [at this website](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html) or, alternatively, by running:
```
 python ./bin/unique-kmers.py -k <your-read-length> <path/to/genome/fasta/file.fa>
``` 
This script will return you the total estimated number of k-mers found in your species genome assembly. If you need further info on this script look at [MR Crusoe *et al.*, 2015](http://dx.doi.org/10.12688/f1000research.6924.1). <br/>

   - A directory containing a Bowtie2 indexed genome. For this you can either run:
```
module load bowtie2/2.4.4
bowtie2-build <your-species-genome>.fa <index-prefix-name>
```
or you could donwload it from the [bowtie webpage](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#:~:text=Bowtie%202%20indexes%20the%20genome,and%20paired%2Dend%20alignment%20modes). In any case, remember to also download from the UCSC, the species genome sequence in the 2bit format as:

```
cd path/to/your/bowtie2/index/genome
wget  http://hgdownload.cse.ucsc.edu/goldenpath/<your-species-assembly>/bigZips/<your-species-assembly>.2bit 
```
and specify this path in relative directives within the `config/snakemake-config.yaml` file.

* Create your own python virtual env and install some modules
As indicated in the [Milton documentation](https://wehieduau.sharepoint.com/sites/rc2/SitePages/Installing-software.aspx?Mode=Edit#installing-software) is best practice to create your python virtual envirnoment where you can `pip install` all your packages. To do so, you need to run the following:
```
module load python/3.7.0
virtualenv myvenv
. myvenv/bin/activate
```
Afterwards you can install all your python modules. Here I am assuming you have already created a virtual env called `mypyenv`. If not simply run `virtualenv myvenv`. When running this pipeline either interactively or via `sbatch` ([see below]((#running-the-pipeline))), a line in the `config/module.txt` will activate your `mypyenv` python virtual enviroment. If you already have created one yourself then just change this line to whatever your virtual environment is called. <br/>

Regarding the python modules I am using in this pipeline, make sure to have installed in your virtual environment these modules: `pybedtools`, `argparse`, `pandas` and `subprocess` before running this pipeline.

* Install all required R libraries such as `yaml`, `data.table`, `GenomicAlignments`, `ATACseqQC`, `csaw`, `GenomicFeatures`, `TxDb.<your-species>.UCSC.<your-species-assembly>.knownGene`, `UpSetR`. 
All the other libraries I am using here should be pretty standard.

* Once you know how many samples you are preprocessing, replace the following lines within the `utils/r-utils.R`:
```
qualitative_palette = brewer.pal.info[brewer.pal.info$category == 'qual',]
sample_palette = sample(unlist(mapply(brewer.pal, qualitative_palette$maxcolors, rownames(qualitative_palette))),length(samples))
names(sample_palette) = samples
```
with something like:
```
sample_palette =  c('one','palette','per','sample')
```
in order to define a color scheme for your samples.

## Running the pipeline
Read here to understand how to run this pipeline either interactively or using  `sbatch`.
### Interactively
I suggest running interactive sessions to test/develop this pipeline using a subset of your fastq files to see whether everything runs smoothly and whether all desired results are produced. To do so, you then need to run:

```
utils/subsample-files.sh -i <indir> -o <outdir> -n <numbreads>
e.g., ./utils/subsample-files.sh -i /wehisan/general/user_managed/grpu_jchoi_0/received/AnneMarie/annemarie_250822  \
-o /wehisan/general/user_managed/grpu_jchoi_0/projects/davide/atac-pipeline/data/subsampled-files/ \
-n 3000
```
This script will extract the first `n = number of reads` from all the files listed in the `i = input dir`. You can then run the entire pipeline interactively on this subset of reads. However, if you do so, remember to correctly specify the location of your subsampled fastq files in the `config/snakemake-config.yaml` file, which is defined in the `fastqdir` directive. <br/>

Next, set up an interactive session on SLURM. I personally do it using my custom made script `run-interactively.sh`. This script reads all `salloc` arguments (e.g., mem/time etc..) specied in the `config/cluster_config.yaml` file (you can changed them to whatever you might need) and requests those resources without you having to write the entire `salloc` command each time. If you also want to use this utility then move the script outside the project directory and set up your interactive session as:

```
chmod +x run-interactively.sh # it should already be executable though
./run-interactively.sh -p <your-project-name> 
eval "$(cat "$modules")"  # this loads all the modules listed in the config/modules.txt file
```
Now you should have all the resources allocated and modules loaded in your own interactive environment. <br/>

To then run the snakemake pipeline interactively simply type:
```
snakemake --cores 9 -s snakefile-preprocess.smk
```
and once its finished, if it doesnt crash, you should get all your results in the `out/` directory and log files in the `logs/` directory. <br/>

**PS:** Using the subsampled test files the pipeline is pretty quick in finishing. Only the set of rules using deepTools do take quite a bit to complete. This means that with the full fastq files this will likely take a while to finish. However, because 1) the pipeline is modular and 2) I defined a main and a qc group of analyses with different priorities, you will be able to still investigate your peaks of chromatin accessibility while deepTools are still running in the background.

### In the background
To run the pipeline on the background run the following lines:
```
mkdir slurm-report 
sbatch utils/run-on-cluster.sh
```
**PS:**  Be sure to correctly specify the name AND location of your python virtual environment in the corresponding line in the `utils/run-on-cluster.sh` file! There should be no need to change anything else in there.
----
## Pipeline overview
The pipeline contained in this repo goes through the following steps:
1. FastQC
2. Adaptor trimming (Trimmomatic)
3. Alignment (Bowtie2)
4. Post-alignment filtering 
5. Peak-calling (MACS2)
7. QCs generation, e.g. using deepTools and other custom scripts
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
Here I am listing all the QCs that are generated with this pipeline. 
QCs:
* Bowtie2 alignment summary results (e.g., % mapped and uniquely mapped reads)
* Number of reads in each bam file
* Library complexity estimate (a table containing different metrics including PBC1/PBC2 and NRF)
* Number of peaks (following removal of those overlapping blacklisted regions)
* Extra peak QCs (e.g. plot number of peaks by size per sample)
* Fraction Reads in Peaks (FRiP)
* Enrichment of peaks around transcription start sites (TSS)

In addition to these, I am using the [deepTools suite of commands](https://deeptools.readthedocs.io/en/develop/) to further explore the ATAC-seq data and generating other QCs. Here below you can read about some of the results you would get.
#### BAM Summary and Coverage
Both of these QCs are obtained by running deepTools programs. Bam coverage is calculated as the number of reads over short consecutive counting windows of defined size. The snakemake pipeline will run deepTools `bamCoverage` and `plotCoverage` programs to respectively calculate and plot the coverage for each sample. The plotting command will return 2 plots indicating:
1. The frequencies of the observed read coverages per sample
2. The fraction of the genome covered by >= a given number of reads 
To calculate and plot the correlation between read coverages across samples I am using `multiBamSummary` and `plotCorrelation`, respectively. The first program computes the read coverages over the entire genome and/or specific genomic regions for >=2 BAM files, whereas the second one plots a heatmap with the resulting Peason correlation values.
#### GC bias
This QC step is performed to test the assumption of an expected uniform distribution of sequenced reads across the genome, regardless of their base-pair composition. This assumption can fail during library prep, when PCR enriches for GC-rich fragments. Here I am using the deepTools command `computeGCbias` to first calculate the expected GC profile for the genome of the species of interest, as the distribution of GC-content for DNA fragments of a given size. The program will compare the expectations with the observations, returning a plot for each sample showing the number of reads overlapping genomic regions of increasing GC content and the log2(observed/expected) ratio of reads per GC content region. If there is no GC-bias in your samples then the observed and expected GC profiles would be similar, i.e. ratio of observed/expected GC content per fragment length = 1. If you spot a GC-bias, you could run the `correctGCbias` program, which removes reads from regions of too high-coverage and add reads to regions of low-coverage.
#### Cumulative enrichment (BAM fingerprint)
For this QC, I am using the deepTools `plotFingerprint` program. It samples the indexed BAM files and calculates the samples cumulative enrichment by counting all reads overlapping a window (bin) of a given length. The resulting plot is useful to assess how well the signal of the sequenced reads can be differentiated from the background distribution of reads in the control input sample. An ideal input with perfect uniform distribution of reads along the genome (i.e. without enrichments in specific open chromatin regions) and infinite sequencing coverage would generate a straight diagonal line. On the other hand, a strong enrichment in very specific (i.e. fragment size) chromatin regions will result in a prominent and steep rise of the cumulative sum towards the highest rank. This means that a big chunk of reads from the test (ChIP/ATAC) sample is located in few bins which corresponds to high, narrow enrichments typically seen for transcription factors. <br/>

**PS**: check where cumulative curve starts on the plots as this will give you an indication of the percentage of the genome that was not sequenced at all (i.e. bins containing 0 reads). However, such percentage can be quite high, especially for extremely high enrichment which would result in the vast majority of reads occurring within few peaks, [see these examples](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html).
