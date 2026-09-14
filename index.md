# Metabarocding practical - using soil samples from Silwood Park

## Introduction

In this practical, we will go from raw sequences to community-level analyses. The samples used in this study were collected from the [Ecological Fractal Network](https://ecofracnetwork.github.io) points at Silwood Park. You can see the specific collection points [here](https://www.google.com/maps/d/viewer?mid=1gYaoOn5ypAK2B-bL8uXjdSTMu4CWCyI&ll=51.40958062885105%2C-0.6467416439178342&z=15).

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

Here we are using **16S rRNA gene amplicon sequencing data generated on an Illumina MiSeq platform, targeting the bacterial community**.

### Using R

This practical relies on a basic understanding of the R language and R studio. I recommend that you create a new RStudio environment to run this analysis within.

## Task 1: Download data

```r

# Set your path - this is where your data will be downloaded to

path <- "path/to/somewhere/on/your/computer/or/OneDrive"

# Download data

wget https://raw.githubusercontent.com/theobrook/Silwood_Park_Soil_Metabarcoding_Practical/main/data/reads.fastq.gz

# Check the data downloaded successfully

list.files(path)

# Set your save_path, this is where all your outputs will be saved

save_path <- "path/to/somewhere/on/your/computer/or/OneDrive" # Not the same place as your path above

```

## Task 2: Install and load packages (libraries)

```r
# Install libraries

install.packages(data2)
install.packages(ggplot2)

# Load libraries

library(dada2)
library(ggplot2)

```

## Task 3: Load data

```r
# Forward and reverse fastq filenames have format: SAMPLENAME_1.fq and SAMPLENAME_2.fq

fnFs <- sort(list.files(path, pattern="_1.fq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fq", full.names = TRUE))

# Extract sample names

sample.names <- sapply(strsplit(basename(fnFs), "_"), function(x) paste(x[-length(x)], collapse = "_"))

```

## Task 4: Inspect read quality profiles ##

```r
# Forward reads

quality_profiles_fnFs <- plotQualityProfile(fnFs[1:66])

# Save - each sample as a separate PNG

for (i in 1:31) {
  p <- plotQualityProfile(fnFs[i])
  ggsave(filename = file.path(save_path, paste0("quality_profiles/quality_profile_forward_", i, ".png")), 
         plot = p, width = 10, height = 7)
}

# Reverse reads
quality_profiles_fnRs <- plotQualityProfile(fnRs[1:66])

# Save - each sample as a separate PNG
for (i in 1:31) {
  p <- plotQualityProfile(fnRs[i])
  ggsave(filename = file.path(save_path, paste0("quality_profiles/quality_profile_reverse_", i, ".png")), 
         plot = p, width = 10, height = 7)
}

```
