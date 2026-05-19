# =============================================================================
# Script: 04_data_wrangling.R
# Author: Isis Sebastião
# Description: Data wrangling operations in R — binding, handling missing
#              values, categorizing continuous variables, and merging datasets.
#              Examples use simulated genomic and gene expression data.
# =============================================================================

rm(list = ls())
library(dplyr)

# =============================================================================
# PART 1: Binding and correlation
# =============================================================================

# Simulated expression data: two genes across samples
gene_A_expr <- c(5.2, 8.1, 12.4, 3.3, 9.0, 6.7, 11.2, 4.5)
gene_B_expr <- c(4.9, 7.8, 13.1, 3.0, 8.7, 6.2, 10.9, 4.1)

# Combine into a data frame using cbind
expr_matrix <- as.data.frame(cbind(gene_A = gene_A_expr,
                                   gene_B = gene_B_expr))
head(expr_matrix)

# Pearson correlation between two genes
cor(expr_matrix$gene_A, expr_matrix$gene_B)
# A high positive correlation may indicate co-expression or co-regulation.

# Logical check: verify data integrity after binding
gene_A_expr == expr_matrix$gene_A

# Sort expression values (ascending)
gene_A_sorted <- sort(gene_A_expr)

# =============================================================================
# PART 2: Handling missing data (NA values)
# =============================================================================

# Simulated metadata table with missing values
metadata <- data.frame(
  sample_id  = paste0("S", 1:8),
  condition  = c("control", "treatment", "control", "treatment",
                 "control", "treatment", "control", "treatment"),
  read_depth = c(25e6, 30e6, NA, 28e6, 22e6, NA, 31e6, 27e6),
  rin_score  = c(8.5, 7.9, 8.1, NA, 7.5, 8.0, NA, 7.8)
)

# Detect missing values
is.na(metadata)
colSums(is.na(metadata))   # count NAs per column

# Replace a known erroneous value with NA (e.g., sentinel value 999)
metadata$read_depth[3] <- 999
metadata$read_depth     <- na_if(metadata$read_depth, 999)

# Remove rows with any NA (use carefully — may reduce sample size)
metadata_clean <- na.omit(metadata)
dim(metadata_clean)

# Calculate mean ignoring NAs
mean(metadata$read_depth, na.rm = TRUE)
mean(metadata$rin_score,  na.rm = TRUE)

# =============================================================================
# PART 3: Categorizing continuous variables
# =============================================================================

# Simulated log2 fold-change values for a gene set
set.seed(42)
log2fc <- round(rnorm(100, mean = 0, sd = 2), 3)

# Quartile thresholds
quantile(log2fc, na.rm = TRUE)
quantile(log2fc, 0.75, na.rm = TRUE)

# Categorize using cut() — 4 expression groups based on quartiles
q <- quantile(log2fc, na.rm = TRUE)
fc_category <- cut(
  log2fc,
  breaks = c(-Inf, q[2], q[3], q[4], +Inf),
  labels = c("Very Low", "Low", "High", "Very High")
)
table(fc_category)

# Categorize using ifelse() — biologically meaningful thresholds
# Standard cutoffs: |log2FC| > 1 = DE, padj < 0.05 = significant
fc_direction <- ifelse(log2fc > 1,  "upregulated",
                ifelse(log2fc < -1, "downregulated",
                                    "not significant"))
table(fc_direction)

# Convert to factor with explicit level order
fc_direction <- factor(fc_direction,
                       levels = c("downregulated", "not significant", "upregulated"))
table(fc_direction)

# =============================================================================
# PART 4: Merging datasets
# =============================================================================

# Gene expression results (e.g., DESeq2 output)
de_results <- data.frame(
  gene_id = c("AT1G01010", "AT1G01020", "AT1G01030", "AT1G01040"),
  log2fc  = c(2.1, -1.4, 0.3, 3.8),
  padj    = c(0.001, 0.043, 0.890, 0.0001)
)

# Functional annotation table (e.g., from InterProScan or MAKER)
annotation <- data.frame(
  gene_id    = c("AT1G01010", "AT1G01020", "AT1G01030", "AT1G01050"),
  go_term    = c("GO:0006950", "GO:0003700", "GO:0005515", "GO:0006355"),
  description = c("stress response", "transcription factor",
                  "protein binding",  "gene expression regulation")
)

# Merge by gene_id (inner join: keeps only matching rows)
merged_inner <- merge(de_results, annotation, by = "gene_id")
head(merged_inner)

# Merge keeping all DE genes (left join)
merged_left <- merge(de_results, annotation, by = "gene_id", all.x = TRUE)
head(merged_left)

# Merge with different column names in each table
annotation$locus_id <- annotation$gene_id
merged_renamed <- merge(de_results, annotation,
                        by.x = "gene_id", by.y = "locus_id")
head(merged_renamed)

# =============================================================================
# PART 5: Transposing and reshaping
# =============================================================================

# Expression matrix: genes as rows, samples as columns (standard format)
expr_wide <- data.frame(
  gene_id  = c("AT1G01010", "AT1G01020", "AT1G01030"),
  ctrl_1   = c(100, 250, 80),
  ctrl_2   = c(110, 230, 90),
  treat_1  = c(300, 120, 85),
  treat_2  = c(280, 130, 78)
)

# Transpose: samples as rows, genes as columns
expr_t <- as.data.frame(t(expr_wide[, -1]))  # exclude gene_id column
colnames(expr_t) <- expr_wide$gene_id
expr_t <- tibble::rownames_to_column(expr_t, var = "sample_id")
print(expr_t)

# Rename a column
colnames(expr_wide)[colnames(expr_wide) == "ctrl_1"] <- "control_replicate_1"
