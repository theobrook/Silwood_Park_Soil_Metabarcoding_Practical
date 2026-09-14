# Silwood Park Soil Metabarcoding Practical

## Introduction

In this practical, we will go from raw sequences to community-level analyses. The samples used in this study were collected from the [Ecological Fractal Network](https://ecofracnetwork.github.io) points at Silwood Park. You can see the specific collection points [here](https://www.google.com/maps/d/viewer?mid=1gYaoOn5ypAK2B-bL8uXjdSTMu4CWCyI&ll=51.40958062885105%2C-0.6467416439178342&z=15).

This practical is an adaptation of the [DADA2 Pipeline Tutorial (1.16)](https://benjjneb.github.io/dada2/tutorial.html), a useful resource for beginners in microbial bioinformatics.

The specific aims are to:
1. Understand the output of a short-read sequencing machine
2. Perform quality control on raw sequencing reads
3. Perform community-level analyses using amplicon sequence variants (ASVs)
4. Assign taxonomy to ASVs, creating virtual taxa (VTs)
5. Love bioinformatics!

### Key terms

- **Amplicon sequencing (metabarcoding)** = a technique where a specific, short region of DNA (the "marker gene") shared by many organisms is amplified via PCR and sequenced from a mixed environmental sample. This lets us profile every organism present at once rather than sequencing one species at a time.
- **Short-read sequencing** = a sequencing technology (e.g., Illumina platforms) that reads DNA in short fragments (usually 50–300 base pairs long). Although "long-read" technologies which read thousands of bases at once are available, short reads are cheaper but require more computational work to reconstruct the full picture.
- **Raw sequencing reads** = the unprocessed output straight off the sequencing machine, usually stored as FASTQ files, containing both the DNA sequence and a quality score for each base call.
- **Quality control (QC)** = the process of filtering and trimming raw reads to remove low-quality base calls, sequencing adapters, and other artefacts before any biological analysis, so downstream results reflect real biology rather than sequencing noise.
- **Amplicon sequence variant (ASV)** = a unique DNA sequence recovered from your samples, inferred at single-nucleotide resolution. ASVs are the modern, higher-resolution replacement for the older approach of clustering reads into OTUs (operational taxonomic units).
- **Taxonomy assignment** = matching each ASV against a reference database to determine which organism (or group of organisms) it most likely came from.
- **Virtual taxa (VT)** = clusters of ASVs that represent the same underlying taxon once taxonomy has been assigned, used to consolidate sequence-level variation into biologically meaningful groups for community analysis.

Here we are using 16S rRNA **TBC!!!** gene amplicon sequencing data generated on an Illumina MiSeq platform, targeting the bacterial community.

### Using R

This practical relies on a basic understanding of the R and RStudio. I recommend that you create a new RStudio environment to run this analysis within.

## Task 1: Download data

We extracted DNA from the soil samples using the [Qiagen DNeasy PowerSoil Pro Kit](https://www.qiagen.com/us/products/discovery-and-translational-research/dna-rna-purification/dna-purification/microbial-dna/dneasy-powersoil-pro-kit) and sequenced the 16S rRNA gene. **TBC!!!** There is a fair amount of data, so you might want to consider downloading it to your Imperial OneDrive account.

```r
# Set your path - this is where your data will be downloaded to
path <- "path/to/somewhere/on/your/computer/or/OneDrive"

# Download data
wget https://raw.githubusercontent.com/theobrook/Silwood_Park_Soil_Metabarcoding_Practical/main/data/reads.fastq.gz # if you struggle with wget, you can download it manually to your desired folder

# Check the data downloaded successfully
list.files(path)

# Set your save_path, this is where all your outputs will be saved
save_path <- "path/to/somewhere/on/your/computer/or/OneDrive" # do not use the same location as your path above
```

## Task 2: Install and load packages (libraries)

```r
# Install libraries
install.packages("data2") # a bioinformatic package to denoise amplicon sequencing data and infer ASVs
install.packages("phyloseq") # a bioinformatic package to import, store, analyse, and plot microbiome (and phylogenetic) sequencing data
install.packages("ggplot2") # a package for plotting

# Load libraries
library(dada2)
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

If `length(fnFs)` and `length(fnRs)` don't match, it usually means a forward or reverse file is missing for one sample. Double check your `data` folder before continuing.

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

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(save_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(save_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Filter and trim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(240, 210),
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
sample.names <- sapply(strsplit(basename(filtFs), "_"), 
                       function(x) paste(x[1:(length(x)-3)], collapse="_"))

