IDEAFIX
=======

ideafix is a decision tree-based variant refinement tool that filters
formaldehyde-induced cytosine deaminations from variant lists obtained
from DNA sequencing data from FFPE specimens.

Installation
------------

``` r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("mmaitenat/ideafix")
```

Requirements
------------

ideafix needs the following programs to run:

-   bcftools

-   samtools

ideafix also needs the following files: \* Fasta file of the genome data
was aligned to. If data is unaligned, the genome to align the data to.

If variant calling is to be run from bam files, you will also need:

-   GATK4

Example
-------

``` r
library(ideafix)
```
