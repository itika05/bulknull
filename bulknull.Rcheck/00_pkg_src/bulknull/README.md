# bulknull: Post-Hoc Diagnostic for Dilution in Cell-Type-Specific Expression

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/dritikaarora/bulknull/graph/badge.svg)](https://app.codecov.io/gh/dritikaarora/bulknull)
<!-- badges: end -->

An R package for diagnosing whether an observed bulk-level null result reflects true dilution of a single-cell signal (small cluster fraction + hypothesised within-cluster effect) versus a false-positive single-cell call.

**Important:** bulknull asks "IF there were a within-cluster IPF-vs-control effect of size X, would we have power to detect it?" It does NOT measure within-cluster effects directly. You supply a hypothesised effect size; the package computes power and posterior probability of dilution.

## Overview

`bulknull` provides two diagnostic quantities:

1. **DDS (Dilution Discordance Score)**: Bayesian posterior probability that the observed bulk-level null reflects true dilution (values 0–1, ~0.5 is ambiguous)

2. **DI (Diagnosability Index)**: Statistical power to detect dilution if it truly occurred at the observed single-cell effect size (values 0–1, <0.3 is underpowered)

## Installation

```r
# From GitHub (when available)
# devtools::install_github("arora-lab/bulknull")

# From local source
devtools::load_all("path/to/bulknull")
```

## Quick Start

**Primary use case: Study design (required_design)**

If you want to know what study design (cluster fraction or bulk measurement precision) is needed to achieve a target statistical power:

```r
library(bulknull)

# Example: What cluster fraction is needed to achieve DI = 0.5 (50% power)?
design <- required_design(
  target_di = 0.5,
  sc_zscore = 2.077409,      # Module enrichment score (M8 vs. other myeloid)
  condition_fraction = 0.631,  # Case fraction in cohort
  n_cluster = 1332,            # Cells in target cluster
  cluster_fraction = 0.069,    # Current cluster fraction
  bulk_se = 0.062774           # Current bulk measurement SE
)
# Result: requires cluster_fraction = 0.8753 (impossible) or bulk_se = 0.0049818 (159x larger cohort)
```

**Secondary use case: Post-hoc diagnosis (bulknull)**

Given an observed bulk null and single-cell signal, diagnose whether dilution is likely:

```r
# Example: DPP9 in M8 cells (k=5 tissue-only)
result <- bulknull(
  bulk_beta = -0.084210,
  bulk_se = 0.062774,
  bulk_fdr = 0.449745,
  sc_zscore = 2.077409,
  cluster_fraction = 0.069,
  condition_fraction = 0.631,
  n_cluster = 1332
)

# Interpretation:
# - DDS = 0.454 (≈50% posterior probability of dilution; ambiguous)
# - DI = 0.065 (6.5% power; severely underpowered to detect dilution)
# - Verdict: UNDERPOWERED (study lacks power to distinguish dilution from true null)
```

## Obtaining sc_zscore (Hypothesised Within-Cluster Effect)

The sc_zscore input is NOT a measured IPF-vs-control effect within the cluster. Instead, it is a hypothesised magnitude for such an effect, used to ask study design questions.

**To obtain sc_zscore, compute a donor-level pseudobulk IPF-vs-control z-score within your target cluster:**

1. Subset single-cell expression to your target cluster (e.g., M8 myeloid cells only)
2. For each donor, compute mean expression across cells in that cluster (pseudobulk)
3. Perform donor-level differential expression: IPF case donors vs. control donors
4. Extract the z-score (or t-statistic converted to z) for your gene of interest
5. Use this donor-level z-score as sc_zscore input

**Do NOT use:**
- Between-state enrichment scores (e.g., M8 vs. other myeloid clusters). These measure cell-type distinctiveness, not case-vs-control effects.
- Cell-level Welch t-tests. These confound within-cluster heterogeneity with case-vs-control signal.

**Example (pseudocode):**
```r
# Subset to target cluster
m8_cells <- seurat_obj[, seurat_obj$cluster == "M8"]

# Pseudobulk: mean expression per donor
pseudobulk <- aggregate_by_donor(m8_cells)  # rows=genes, cols=donors, values=mean expr

# Donor-level DE (case vs. control)
case_donors <- pseudobulk[, metadata$condition == "case"]
ctrl_donors <- pseudobulk[, metadata$condition == "control"]

# T-test and convert to z-score
t_stat <- t.test(case_donors["GENEX", ], ctrl_donors["GENEX", ])$statistic
sc_zscore <- sign(t_stat) * sqrt(t_stat^2 / df)
```

## Framework Requirements

The bulknull framework is only applicable when:
- **Bulk-level result is null** (FDR > 0.05 or specified threshold)
- **Target cluster shows strong signal** in some context (e.g., enriched vs. other clusters, or case-specific effect)
- **Cluster fraction is small** (typically 0.01–0.25)

If the bulk effect is significant (FDR < 0.05), use standard cell-type-specific DE methods instead (CARseq, TOAST, bMIND, etc.).

## How It Works

### DDS (Dilution Discordance Score)

DDS computes the Bayesian posterior probability that the observed bulk z-score came from the dilution hypothesis (H1: z ~ N(μ, 1)) versus the null (H0: z ~ N(0, 1)):

```
DDS = φ(z_bulk − μ) / (φ(z_bulk − μ) + φ(z_bulk))
```

where φ is the standard normal PDF and μ_dilution = (cluster_frac × sc_zscore) / (bulk_se × √n_eff).

Interpretation:
- **DDS > 0.7**: Dilution probable
- **DDS ≈ 0.5**: Ambiguous
- **DDS < 0.3**: Null probable (signal is truly absent in bulk)

### DI (Diagnosability Index)

DI estimates the probability of correctly detecting dilution if it occurred:

```
DI = Φ(μ − z_critical)
```

where z_critical = qnorm(0.95) = 1.644854 for one-sided α=0.05, and Φ is the standard normal CDF.

Power categories:
- **DI < 0.3**: Low power (undetectable)
- **DI 0.3–0.8**: Moderate power (unreliable)
- **DI > 0.8**: High power (reliable detection)

## Important Limitation

At realistic cluster fractions (0.05–0.15), DI did not exceed 0.10 in our portfolio (one-sided α=0.05), making dilution claims **effectively unfalsifiable** with RNA-seq data alone. This motivates spatial profiling (e.g., Xenium, spatial proteomics) as validation.

## Testing

The package includes 129 automated test assertions covering function correctness, input validation, applicability gating, formula invariants, and validation against independent implementations. Documentation-code consistency for critical values is verified manually; see NEWS.md.

**Before committing manuscript changes:** Run `Rscript scripts/check_orphan_numbers.R` to ensure all numeric values in `.md` files appear in `NUMBERS_TABLE.tsv`. This guard test inspects project markdown for orphan numbers and rejects any value not in the approved list.

## Citation

If you use bulknull, please cite:
> Arora I, Yaqinuddin A. (2026). bulknull: Post-hoc diagnostic for dilution in cell-type-specific expression. *bioRxiv* (preprint) / *Bioinformatics* (if published).

## License

MIT License. See LICENSE file.

## Authors

- Itika Arora
- Alfaisal University

## References

- bulknull Application Note (in preparation)
- Phase 0 synthetic validation & Phase 1 consistency check (included in vignette)
