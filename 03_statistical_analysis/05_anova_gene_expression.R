# =============================================================================
# Script: 05_anova_gene_expression.R
# Author: Isis Sebastião
# Description: One-way ANOVA and non-parametric alternative (Kruskal-Wallis)
#              applied to relative gene expression data (e.g., from qRT-PCR).
#              Includes assumption testing (normality, homogeneity of variance)
#              and post-hoc analysis (Tukey HSD).
#
# Biological context: comparing relative expression of 6 candidate genes
#              across treatment groups — a common step in validating
#              RNA-seq differential expression results.
# =============================================================================

# --- 0. Dependencies ---------------------------------------------------------

# install.packages("car")
# install.packages("ggplot2")
library(car)
library(ggplot2)

# --- 1. Input data -----------------------------------------------------------

# Relative expression values from qRT-PCR (3-4 biological replicates per gene)
# Values represent 2^(-deltaCt) normalized to reference gene
genes_expr <- data.frame(
  gene = rep(c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"), each = 4),
  expression = c(
    2.7, 0.4, 0.6, 1.4,   # GeneA
    1.6, 1.1, 0.7, 0.8,   # GeneB
    1.9, 1.5, 0.5, 1.0,   # GeneC
    1.7, 1.6, 0.5, 0.6,   # GeneD
    0.7, 1.1, 1.0, 1.0,   # GeneE
    2.4, 1.0, 0.8, 0.4    # GeneF
  )
)

# Convert gene to factor (required for ANOVA)
genes_expr$gene <- factor(genes_expr$gene)

# --- 2. Exploratory visualization --------------------------------------------

ggplot(genes_expr, aes(x = gene, y = expression, fill = gene)) +
  geom_boxplot(alpha = 0.7, color = "black") +
  geom_jitter(width = 0.15, alpha = 0.8, size = 2) +
  labs(
    title = "Relative expression of candidate genes",
    subtitle = "Values normalized to reference gene (2^-deltaCt)",
    x = "Gene",
    y = "Relative expression"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# --- 3. One-way ANOVA --------------------------------------------------------

anova_result <- aov(expression ~ gene, data = genes_expr)
summary(anova_result)

# Interpretation:
# If p < 0.05, at least one gene differs significantly from the others.
# The ANOVA does NOT tell you which genes differ — use Tukey HSD for that.

# --- 4. Assumption testing ---------------------------------------------------

# 4a. Normality of residuals (Shapiro-Wilk test)
residuals_anova <- residuals(anova_result)
shapiro_result  <- shapiro.test(residuals_anova)
print(shapiro_result)
# H0: residuals are normally distributed
# If p > 0.05, normality is not rejected

# Q-Q plot to visually assess normality
qqnorm(residuals_anova, main = "Q-Q plot of residuals")
qqline(residuals_anova, col = "red", lwd = 2)

# 4b. Homogeneity of variances (Levene's test)
levene_result <- leveneTest(expression ~ gene, data = genes_expr)
print(levene_result)
# H0: variances are equal across groups
# If p > 0.05, homogeneity of variance is not rejected

# --- 5. Decision: ANOVA or Kruskal-Wallis ------------------------------------

normality_ok <- shapiro_result$p.value > 0.05
homogeneity_ok <- levene_result$`Pr(>F)`[1] > 0.05

if (normality_ok & homogeneity_ok) {

  cat("\n✅ Assumptions met: using parametric ANOVA\n\n")
  print(summary(anova_result))

  # Post-hoc: Tukey HSD — identifies which gene pairs differ significantly
  tukey_result <- TukeyHSD(anova_result)
  print(tukey_result)
  plot(tukey_result, las = 1)

} else {

  cat("\n⚠️ Assumptions not met: using non-parametric Kruskal-Wallis\n\n")
  kruskal_result <- kruskal.test(expression ~ gene, data = genes_expr)
  print(kruskal_result)

  # Post-hoc for Kruskal-Wallis: pairwise Wilcoxon with correction
  pairwise.wilcox.test(genes_expr$expression, genes_expr$gene,
                       p.adjust.method = "BH")
}

# --- 6. Notes ----------------------------------------------------------------

# This workflow is analogous to comparing expression across accessions
# in pangenome studies — where expression variation between genotypes
# needs to be statistically assessed before being attributed to
# structural genomic differences (e.g., PAV — presence/absence variation).
