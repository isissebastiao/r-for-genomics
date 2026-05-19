# =============================================================================
# Script: 02_vectors_and_data_types.R
# Author: Isis Sebastião
# Description: R data types and vector operations illustrated with
#              genomics-relevant examples (gene names, expression values,
#              chromosomal positions).
# =============================================================================

rm(list = ls())

# --- 1. Numeric vectors ------------------------------------------------------

# Gene expression values (log2 fold-change)
log2fc <- c(2.1, -1.5, 0.8, 3.4, -0.3)
print(log2fc)

# Chromosomal positions of SNPs
snp_positions <- c(1000, 2500, 4800, 10200, 30500)

# Generating sequences (e.g., genomic coordinates at regular intervals)
coord_window <- seq(from = 1000, to = 2000, by = 100)

# Repeating values (e.g., replicates per treatment)
replicates <- rep(1:3, times = 4)     # 4 treatments, 3 replicates each
replicates_each <- rep(1:3, each = 4) # each replicate ID repeated 4x

# --- 2. Vector operations ----------------------------------------------------

# Sort expression values (ascending)
sort(log2fc)

# Reverse order
rev(log2fc)

# Frequency table (e.g., SNP counts per chromosome)
chromosomes <- c("Chr1", "Chr2", "Chr1", "Chr3", "Chr1", "Chr2")
table(chromosomes)

# Unique values
unique(chromosomes)

# Subsetting vectors
log2fc[1]              # first element
log2fc[-1]             # all except first
log2fc[1:3]            # elements 1 to 3

# Logical subsetting — selecting differentially expressed genes (|log2FC| > 1)
log2fc[abs(log2fc) > 1]

# Using %in% — check which chromosomes are in a set of interest
chromosomes[chromosomes %in% c("Chr1", "Chr3")]

# --- 3. Character vectors ----------------------------------------------------

# Gene identifiers
gene_ids <- c("AT1G01010", "AT1G01020", "AT1G01030", "AT2G00010")

# Sample names
samples <- c("control_1", "control_2", "treatment_1", "treatment_2")

# --- 4. Factor vectors -------------------------------------------------------

# Treatment groups as factors (important for statistical models)
treatment <- factor(
  c("control", "treatment", "control", "treatment"),
  levels = c("control", "treatment")   # reference level = control
)
print(treatment)

# --- 5. Data type checking and conversion ------------------------------------

# Check the type of a vector
class(log2fc)       # "numeric"
class(gene_ids)     # "character"
class(treatment)    # "factor"
mode(log2fc)        # "numeric"

# Convert character vector of numbers to numeric
count_values_chr <- c("100", "250", "80", "300")
count_values_num <- as.numeric(count_values_chr)
print(count_values_num)

# --- 6. Mathematical functions -----------------------------------------------

expression_values <- c(5.2, 8.7, 2.1, 15.4, 0.3)

log(expression_values)        # natural log
log(expression_values, 2)     # log2 (standard in transcriptomics)
log(expression_values, 10)    # log10

exp(1)                        # inverse of natural log

round(2.718281, digits = 3)   # rounding

mean(expression_values)
median(expression_values)
sd(expression_values)         # standard deviation
var(expression_values)        # variance
quantile(expression_values)   # quartiles
