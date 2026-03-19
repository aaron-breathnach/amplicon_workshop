# Pre-workshop setup

Please complete these steps before the workshop.

## 1. Install R and RStudio

Install a recent version of:

- R
- RStudio Desktop

## 2. Install required packages

Run the following in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("dada2", "phyloseq", "Biostrings"))
install.packages(c("ggplot2", "dplyr", "readr", "tibble"))
```

## 3. Test your installation

Run:

```r
library(dada2)
library(phyloseq)
```

If both packages load without errors, your setup should be ready.

## 4. Download workshop materials

Clone or download this repository and open the `scripts/amplicon_dada2_workflow.R` script in RStudio.
