# IDEAFIX

ideafix is a decision tree-based variant refinement tool that filters
formaldehyde-induced cytosine deaminations from variant lists obtained
from DNA sequencing data from FFPE specimens.

## Installation

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mmaitenat/ideafix", build_vignettes = TRUE)
```

## Requirements

ideafix needs the following programs to run:

  - bcftools

  - samtools

ideafix also needs the following files:

  - VCF file obtained from an FFPE specimen. Variant calling needs to be
    run with Mutect2, with strand bias annotation enabled. It can be
    either a tumor-only or a normal-tumor paired variant calling.

  - Fasta file of the genome data was aligned to. If data is unaligned,
    the genome to align the data to.

If variant calling is to be run from bam files (feature not available
yet), you will also need:

  - GATK4

## Example

``` r
library(ideafix)
# Extract C:G > T:A variants with a VAF < 30% and their descriptors from vcf filename
vcf_filename <- "mysample_chr15.vcf"
ref_genome <- "chr15.fasta"
descriptors <- get_descriptors(vcf_filename = vcf_filename, fasta_filename = ref_genome)
# Filter variants
predictions_XGBoost <- classify_variants(variant_descriptors = descriptors, algorithm = "XGBoost")
```
