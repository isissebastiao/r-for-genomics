# =============================================================================
# Script: 03_dataframes.R
# Author: Isis Sebastião
# Description: Data frame operations in R using a simulated gene annotation
#              table as example — reflecting common data structures in
#              genomics workflows (e.g., output from functional annotation tools).
# =============================================================================

rm(list = ls())

# --- 1. Reading tabular data -------------------------------------------------

# Reading a tab-separated file (common format for bioinformatics outputs)
# df <- read.table("data/gene_annotation.txt", header = TRUE, sep = "\t")

# Reading a CSV
# df <- read.csv("data/expression_matrix.csv", header = TRUE)

# Reading an Excel file
# library(readxl)
# df <- read_excel("data/samples_metadata.xlsx", sheet = "metadata")

# --- 2. Creating a data frame from scratch (example) -------------------------

# Simulated gene annotation table (typical output from MAKER or BRAKER)
df <- data.frame(
  gene_id    = c("AT1G01010", "AT1G01020", "AT1G01030", "AT1G01040", "AT1G01050"),
  chr        = c("Chr1", "Chr1", "Chr1", "Chr1", "Chr1"),
  start      = c(3631, 5928, 6788, 11649, 23121),
  end        = c(5899, 8737, 9130, 13714, 31227),
  strand     = c("+", "-", "+", "+", "-"),
  log2fc     = c(2.1, -1.4, 0.3, 3.8, -2.6),
  padj       = c(0.001, 0.043, 0.890, 0.0001, 0.012),
  annotation = c("MYB transcription factor", "Stress response protein",
                 "Unknown", "LTR retrotransposon", "Kinase")
)

# --- 3. Exploring a data frame -----------------------------------------------

dim(df)          # dimensions: rows x columns
names(df)        # column names (same as colnames)
colnames(df)
class(df)

head(df)         # first 6 rows
tail(df)         # last 6 rows

# --- 4. Subsetting -----------------------------------------------------------

# By position: [rows, columns]
df[1:3, 1:4]

# By column name
df$log2fc
df[, "annotation"]

# Specific rows and columns
df[c(1, 3, 5), c("gene_id", "log2fc", "padj")]

# --- 5. Filtering with subset() ----------------------------------------------

# Genes with significant differential expression (padj < 0.05)
sig_genes <- subset(df, padj < 0.05)
dim(sig_genes)

# Upregulated significant genes (log2FC > 1 AND padj < 0.05)
upregulated <- subset(df, log2fc > 1 & padj < 0.05)
print(upregulated[, c("gene_id", "log2fc", "padj", "annotation")])

# Genes located on the positive strand
positive_strand <- subset(df, strand %in% c("+"))
dim(positive_strand)

# --- 6. Summary statistics ---------------------------------------------------

mean(df$log2fc)
table(df$strand)
table(df$chr, df$strand)

# --- 7. Transposing a data frame ---------------------------------------------

df_t <- as.data.frame(t(df))
class(df_t)

# --- 8. Adding and modifying columns -----------------------------------------

# Calculate gene length from coordinates
df$gene_length <- df$end - df$start
head(df[, c("gene_id", "start", "end", "gene_length")])

# Classify expression direction
df$direction <- ifelse(df$log2fc > 0, "upregulated", "downregulated")
table(df$direction)

# --- 9. Column and row binding -----------------------------------------------

# cbind: add columns (must have same number of rows)
rpkm_values <- c(45.2, 12.1, 0.8, 102.4, 3.3)
df_extended <- cbind(df, rpkm = rpkm_values)
head(df_extended)

# rbind: add rows (must have same columns)
new_gene <- data.frame(
  gene_id = "AT1G01060", chr = "Chr1", start = 31500, end = 33000,
  strand = "+", log2fc = 1.2, padj = 0.03,
  annotation = "Heat shock protein", gene_length = 1500,
  direction = "upregulated"
)
df_full <- rbind(df[, names(new_gene)], new_gene)
dim(df_full)
