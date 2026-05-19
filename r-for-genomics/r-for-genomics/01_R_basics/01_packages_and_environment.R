# =============================================================================
# Script: 01_packages_and_environment.R
# Author: Isis Sebastião
# Description: Introduction to R package management and environment setup,
#              with focus on packages commonly used in genomic data analysis.
# =============================================================================

# --- 1. Installing and loading CRAN packages ---------------------------------

# General-purpose packages
install.packages("dplyr")    # data manipulation
install.packages("ggplot2")  # data visualization
install.packages("readxl")   # reading Excel files
install.packages("tidyr")    # data tidying

library(dplyr)
library(ggplot2)

# --- 2. Bioconductor packages (used in genomics/transcriptomics) -------------

# Bioconductor requires its own package manager
install.packages("BiocManager")

# Key packages for genomic analyses
BiocManager::install("DESeq2")           # differential expression analysis
BiocManager::install("edgeR")            # RNA-seq differential expression
BiocManager::install("GenomicRanges")    # genomic interval operations
BiocManager::install("Biostrings")       # biological sequence handling
BiocManager::install("VariantAnnotation") # VCF file handling

library(DESeq2)

# --- 3. Getting help ---------------------------------------------------------

# Access documentation for a loaded package
?ggplot2

# Search all available documentation (broader search)
??DESeq2

# View a function's usage example
example(sum)

# --- 4. Environment management -----------------------------------------------

# Assign values to objects
sample_size <- 5
gene_count  <- 20000

# List all objects in the current environment
ls()

# Remove specific objects to free memory
rm(sample_size)

# Clear the entire environment (use carefully)
rm(list = ls())
