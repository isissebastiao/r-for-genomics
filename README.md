# R for Genomic Data Analysis

A collection of R scripts documenting foundational to intermediate workflows
for genomic and transcriptomic data analysis. Developed as part of structured
R training, with examples adapted to biological and genomics contexts.

---

## Motivation

Proficiency in R is a core requirement for modern bioinformatics. These scripts
consolidate practical skills in data manipulation, statistical testing, and
visualization — all of which are directly applied in genomics workflows such as
differential expression analysis, variant filtering, and pangenome data
interpretation.

---

## Repository structure

```
r-for-genomics/
│
├── 01_R_basics/
│   ├── 01_packages_and_environment.R   # CRAN and Bioconductor package management
│   ├── 02_vectors_and_data_types.R     # Data types, vector operations, log2FC examples
│   └── 03_dataframes.R                 # Data frames with simulated gene annotation tables
│
├── 02_data_wrangling/
│   └── 04_data_wrangling.R             # Binding, NA handling, categorization, merging
│
└── 03_statistical_analysis/
    └── 05_anova_gene_expression.R      # ANOVA / Kruskal-Wallis for qRT-PCR data
```

---

## Scripts overview

| Script | Topic | Key functions / packages |
|--------|-------|--------------------------|
| `01_packages_and_environment.R` | Package management | `install.packages`, `BiocManager`, `library` |
| `02_vectors_and_data_types.R` | Vectors & data types | `c`, `seq`, `rep`, `subset`, `log`, `mean`, `sd` |
| `03_dataframes.R` | Data frames | `read.table`, `subset`, `cbind`, `rbind`, `ifelse` |
| `04_data_wrangling.R` | Data manipulation | `dplyr`, `merge`, `na_if`, `na.omit`, `cut`, `t` |
| `05_anova_gene_expression.R` | Statistical testing | `aov`, `TukeyHSD`, `shapiro.test`, `leveneTest`, `kruskal.test`, `ggplot2` |

---

## Biological context

Examples throughout these scripts use data structures common in genomics:

- **Gene annotation tables** — outputs from MAKER, BRAKER, or similar tools
- **Differential expression results** — log2 fold-change and adjusted p-values
- **Sample metadata** — sequencing depth, RIN score, treatment conditions
- **qRT-PCR expression data** — relative expression (2^-deltaCt) per gene

The ANOVA script (`05_anova_gene_expression.R`) is directly applicable to
comparing gene expression across multiple accessions — a common analytical step
in pangenome studies, where expression variation is linked to structural
genomic differences such as presence/absence variation (PAV).

---

## Dependencies

```r
# CRAN
install.packages(c("dplyr", "ggplot2", "readxl", "tidyr", "car", "tibble"))

# Bioconductor
BiocManager::install(c("DESeq2", "GenomicRanges", "Biostrings", "VariantAnnotation"))
```

---

## Author

**Isis Sebastião** — Agronomist | PhD in Biotechnology | Plant Genomics & Bioinformatics  
[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--1596--2523-green)](https://orcid.org/0000-0002-1596-2523)
[![Lattes](https://img.shields.io/badge/Lattes-CNPq-blue)](http://lattes.cnpq.br/5220007563821018)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-isis--sebastiao-blue)](https://www.linkedin.com/in/isis-sebastiao/)
