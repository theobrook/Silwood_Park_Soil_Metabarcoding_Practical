# Bioinformatics using metabarcoding data practical

## Introduction

In this practical, we will go from raw sequences to community-level analyses. The samples used in this study were collected from the [Ecological Fractal Network](https://ecofracnetwork.github.io) points at Silwood Park. You can see the specific collection points [here](https://www.google.com/maps/d/viewer?mid=1gYaoOn5ypAK2B-bL8uXjdSTMu4CWCyI&ll=51.40958062885105%2C-0.6467416439178342&z=15).

The specific aims are to:
1. Understand the output of a **short-read sequencing machine**
2. Perform **quality control** on raw sequencing reads
3. Perform **community-level analyses using amplicon sequence variants (ASVs)**
4. **Assign taxonomy to ASVs**, creating virtual taxa (VTs), and perform basic phylogenetic analysis
5. **Love bioinformatics!**

This practical is an adaptation of the [DADA2 Pipeline Tutorial (1.16)](https://benjjneb.github.io/dada2/tutorial.html), a useful resource for beginners in microbial bioinformatics. If you have any questions, please reach out to Theodore Brook (`T.Brook@kew.org`).

### Key terms

- **Amplicon sequencing (metabarcoding)** = a technique where a specific, short region of DNA (the "marker gene") shared by many organisms is amplified via PCR and sequenced from a mixed environmental sample. This lets us profile every organism present at once rather than sequencing one species at a time.
- **Short-read sequencing** = a sequencing technology (e.g., Illumina platforms) that reads DNA in short fragments (usually 50–300 base pairs long). Although "long-read" technologies which read thousands of bases at once are available, short reads are cheaper but require more computational work to reconstruct the full picture.
- **Raw sequencing reads** = the unprocessed output straight off the sequencing machine, usually stored as FASTQ files, containing both the DNA sequence and a quality score for each base call.
- **Quality control (QC)** = the process of filtering and trimming raw reads to remove low-quality base calls, sequencing adapters, and other artefacts before any biological analysis, so downstream results reflect real biology rather than sequencing noise.
- **Amplicon sequence variant (ASV)** = a unique DNA sequence recovered from your samples, inferred at single-nucleotide resolution. ASVs are the modern, higher-resolution replacement for the older approach of clustering reads into OTUs (operational taxonomic units).
- **Taxonomy assignment** = matching each ASV against a reference database to determine which organism (or group of organisms) it most likely came from.
- **Virtual taxa (VT)** = clusters of ASVs that represent the same underlying taxon once taxonomy has been assigned, used to consolidate sequence-level variation into biologically meaningful groups for community analysis.

## The data

**This practical is currently set up to process one marker type (e.g., 16S, ITS). It can easily be updated if we want to use multiple markers.**

The samples were sequenced on an Illumina MiSeq platform.

### Using R

This practical relies on a basic understanding of the programming language `R`. I recommend that you create a new `RStudio` environment in a `Silwood_Soil_Metabarcoding_Practical” directory to run this analysis within.

### Useful websites
- [Stack Overflow](https://stackoverflow.com/) - a site for programmers to discuss issues/problems
- [DADA2 website](https://benjjneb.github.io/dada2/) - details of the `DADA2` package
- [phyloseq website](https://joey711.github.io/phyloseq/) - details of the `phyloseq` package
- [ggplot2 website](https://ggplot2.tidyverse.org/) - details of the `ggplot2` package

*LLMs such as Claude can be useful, but be careful to check that you understand what they are doing and, crucially, that they are actually doing what you want!*

## Task 1: Download data

We extracted DNA from the soil samples using the [Qiagen DNeasy PowerSoil Pro Kit](https://www.qiagen.com/us/products/discovery-and-translational-research/dna-rna-purification/dna-purification/microbial-dna/dneasy-powersoil-pro-kit) and sequenced the 16S rRNA gene. **TBC!!!** There is a fair amount of data, so you might want to consider downloading it to your Imperial OneDrive account.

The first step is to download the GitHub repository from [here](https://github.com/theobrook/Silwood_Park_Soil_Metabarcoding_Practical) and move the data directory (folder) to somewhere on your device.

Next unzip the directory to extract its contents.

N.B. you might want to use OneDrive as your workspace as there is a large amount of data (**XGB - TBC!!!**).

```r
# Once you have downloaded your data, set your path to where you moved the data directory to. You can find the full file path for a directory by right clicking on the folder and either (a) copying the "Where" field (Mac) or (b) selecting "Properties", and copying the "Location" field (Windows)
path <- "path/to/somewhere/on/your/computer/or/OneDrive"

# (Optional) Tidy up the zip file now that we have extracted its contents
file.remove(zip_dest)

# Check the data downloaded successfully
list.files(path)

# Set your save_path, this is where all your outputs will be saved
save_path <- file.path(pat, "outputs")
```

## Task 2: Install and load packages (libraries)

```r
# Install libraries
install.packages("dada2") # a bioinformatic package to denoise amplicon sequencing data and infer ASVs
install.packages("phyloseq") # a bioinformatic package to import, store, analyse, and plot microbiome (and phylogenetic) sequencing data
install.packages("ggplot2") # a package for plotting

# Load libraries
library(dada2)
library(phyloseq)
library(ggplot2)
```

## Task 3: Load data

Amplicon sequencing such as this usually reads in both directions, creating forward and reverse reads for every DNA fragment. These are stored as two separate files per sample - usually distinguished by a suffix like `_1` (forward) and `_2` (reverse) - and need to be kept paired up, since each forward/reverse pair represents one sequenced fragment.

The data must be read into the R environment. First point R at the folder containing your downloaded FASTQ files, then list and pair up the forward and reverse reads:

```r
# Forward and reverse fastq filenames have format: SAMPLENAME_1.fq and SAMPLENAME_2.fq
fnFs <- sort(list.files(path, pattern="_1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq", full.names = TRUE))

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[-length(x)], collapse = "_"))
```

Before moving on, check that everything has loaded and paired up correctly:

```r
# You should see one name per sample, and the two counts below should match
sample.names
length(fnFs) == length(fnRs)
```

**Checkpoint:** If `length(fnFs)` and `length(fnRs)` don't match, it usually means a forward or reverse file is missing for one sample. Double check your `data` folder before continuing.

## Task 4: Inspect read quality profiles

Before we can filter and trim the reads, we need to know where sequencing quality starts to drop off along each read. DADA2's `plotQualityProfile()` plots this for you.

The grey heatmap shows the frequency of each quality score at each position along the read, the green line is the mean quality score at that position, the orange line is the median, and the orange dashed lines show the 25th and 75th quantiles. As a rule of thumb, quality tends to decline towards the end of the read and you are looking for the position where the mean quality (green line) drops below ~Q30, since that's where you'll want to truncate reads in the next task.

```r
# Set where quality profile plots will be saved
dir.create(file.path(save_path, "quality_profiles"), recursive = TRUE, showWarnings = FALSE)

# Forward reads: quick overview across all samples at once
quality_profiles_fnFs <- plotQualityProfile(fnFs)
quality_profiles_fnFs

# Save each sample's forward-read profile as its own PNG
for (i in seq_along(fnFs)) {
  p <- plotQualityProfile(fnFs[i])
  ggsave(filename = file.path(save_path, "quality_profiles", paste0("quality_profile_forward_", sample.names[i], ".png")),
         plot = p, width = 10, height = 7)
}

# Reverse reads: quick overview across all samples at once
quality_profiles_fnRs <- plotQualityProfile(fnRs)
quality_profiles_fnRs

# Save each sample's reverse-read profile as its own PNG
for (i in seq_along(fnRs)) {
  p <- plotQualityProfile(fnRs[i])
  ggsave(filename = file.path(save_path, "quality_profiles", paste0("quality_profile_reverse_", sample.names[i], ".png")),
         plot = p, width = 10, height = 7)
}
```

**Checkpoint:** Look at your saved quality profiles. At roughly what position do the forward reads start to drop in quality? What about the reverse reads (these are usually a bit worse, can you think of why that might be)? Make a note of these positions, as you'll need them in the next task to set trimming lengths.

## Task 5: Filter and trim reads

Now we have an idea of the quality of our sequences, we need to filter and trim sequences to remove low quality regions.

```r
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(save_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(save_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(X, Y),
                     maxN = 0, 
                     maxEE = c(2, 2), 
                     truncQ = 2, 
                     rm.phix = TRUE,
                     compress = TRUE, 
                     multithread = FALSE)

head(out)
saveRDS(out, file = file.path(save_path, "out.rds"))

# Continue after filtering
filt_path <- file.path(save_path, "filtered")

# List filtered files
filtFs <- sort(list.files(filt_path, pattern="_F_filt.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(filt_path, pattern="_R_filt.fastq.gz", full.names = TRUE))

# Extract sample names from filtered files
sample.names <- sapply(strsplit(basename(filtFs), "_"), function(x) paste(x[1:(length(x)-2)], collapse="_"))
```

**Checkpoint:** You must set `X` and `Y` to the end position you want to truncate the forward and reverse reads to. For example, as these are 250 base pair fragments, to trim the last 30 bases off the forward reads, you would set X to 220.

## Task 6: Learn error rates 

DADA2's core denoising step needs to know what real **sequencing error** looks like before it can distinguish true biological variation and noise. The `learnErrors()` function estimates this by alternating between guessing the error rates and inferring the sample composition until the estimate converges. This producing a model of how likely each type of base-calling error is, at each quality score and at each position along the read. This step can take some time to run, especially on a laptop.

```r
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

# Save the error models so you don't have to re-run this step if you close R
saveRDS(errF, file = file.path(save_path, "errF.rds"))
saveRDS(errR, file = file.path(save_path, "errR.rds"))

# Visualise the estimated error rates
error_plot_F <- plotErrors(errF, nominalQ = TRUE)
ggsave(filename = file.path(save_path, "error_plot_forward.png"), plot = error_plot_F)

error_plot_R <- plotErrors(errR, nominalQ = TRUE)
ggsave(filename = file.path(save_path, "error_plot_reverse.png"), plot = error_plot_R)
```

## Task 7: Sample inference

This is the core denoising step of the DADA2 pipeline. Using the error model learned in **Task 6**, the `dada()` function looks at every unique sequence in each sample and works out which ones represent real biological variants and which are more likely sequencing errors of a more abundant "true" sequence. We now have our **amplicon sequence variants (ASVs)**! 

Note that this is applied separately to the forward and reverse reads. This is usually the slowest step in the whole pipeline, so it will take a while (Windows users especially, since this step doesn't multithread the same way it does on Mac/Linux). **So take a break!**

```r
# Run the core DADA2 denoising algorithm on the forward and reverse reads
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Inspect the returned dada-class object for the first sample
dadaFs[[1]]
dadaRs[[1]]

# Save RDS - this allows you to reload the R data, rather than running the whole script above (e.g. if your R session crashes for whatever reason)
saveRDS(dadaFs, file = file.path(save_path, "dadaFs.rds"))
saveRDS(dadaRs, file = file.path(save_path, "dadaRs.rds"))

# If necessary, you can reload the data if there is an issue
dadaFs <- readRDS(file.path(save_path, "dadaFs.rds"))
dadaRs <- readRDS(file.path(save_path, "dadaRs.rds"))
```

**Checkpoint:** Look at the printed summary for `dadaFs[[1]]` - it should read something like `X sequence variants were inferred from Y input unique sequences.` Why is X usually much smaller than Y? What does that tell you about the relationship between unique sequences and real biological variants?
```

## Task 8: Sample inference
