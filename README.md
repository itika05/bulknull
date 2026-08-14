# bulknull: Bulk-Null Hypothesis Testing for Single-Cell Dilution

An R package implementing Bayesian inference for detecting and quantifying cell-type dilution in bulk tissue transcriptomics. **bulknull** tests whether observed bulk-level effects are driven by true cell-state differences or by compositional dilution—when a cell type of interest represents a small fraction of total tissue.

## The Problem

Bulk tissue RNA-seq measures average gene expression across all cell types. When a cell type of interest is rare (say, 5% of cells), an observed bulk effect could arise from:

1. **True effect**: A large within-cell-type effect d_sc, diluted by the small cluster fraction
2. **Dilution artifact**: Real but small-magnitude signal inflated by compositional effects

**bulknull** quantifies which scenario is more likely, using Bayesian inference and power analysis.

## Key Statistics

### Dilution Discordance Score (DDS)
Posterior probability that observed bulk effect reflects true dilution:
- Range: 0 to 1 (higher = stronger dilution evidence)
- Computed from bulk z-score, cluster proportion, and single-cell effect size

### Diagnosability Index (DI)  
Statistical power to detect dilution at given effect size:
- Range: 0 to 1 (higher = better powered)
- Inverted to solve for required cluster fraction or bulk SE

### Type S and Type M Error Rates
- **Type S**: P(sign error | significant) — risk of wrong-direction conclusion
- **Type M**: E[|z| | sig] / μ — expected magnitude exaggeration when significant

## Installation

```r
# From GitHub
devtools::install_github("itika05/bulknull")
library(bulknull)
```

## Quick Start

```r
# Single gene: DPP9 in M8 lung fibroblasts
result <- bulknull(
  bulk_beta = -0.084,
  bulk_se = 0.063,
  bulk_fdr = 0.45,
  d_sc = 0.118,              # within-cluster effect
  cluster_fraction = 0.07,   # 7% of tissue
  condition_fraction = 0.63, # 63% in condition
  alpha = 0.05
)

print(result)        # Full diagnostics
summary(result)      # One-line summary

# Vectorized: scan 1000+ genes at once
results <- dilution_scan(
  bulk_deg = deg_table,      # data.frame with beta, se, fdr, gene
  sc_zscores = z_vector,     # named numeric vector by gene
  cluster_fraction = 0.07,
  condition_fraction = 0.63,
  n_cluster = 1332           # cells in target cluster
)

head(results)
```

## Core Functions

**Analysis**
- `bulknull()` — Single-gene analysis
- `dilution_scan()` — Vectorized multi-gene analysis
- `dilution_score()` — Compute DDS only

**Design Inversion**
- `required_design()` — What f or bulk_se needed for target DI?
- `critical_precision()` — Minimum bulk SE for target DI?
- `cohort_inflation()` — How much to grow cohort for target DI?
- `dds_bounds()` — Attainable DDS range over z-score span?

**Supporting**
- `diagnosability_index()` — Power given μ, α, sidedness
- `check_dds_applicability()` — Is bulk FDR > threshold?

## Outputs

**bulknull() returns** (list with S3 class):
- `mu_dilution` — Expected z under dilution model
- `w` — Mixture weight: (f·r)/(f·r+(1-f))
- `di` — Diagnosability Index
- `di_ceiling` — Maximum DI when w→1
- `verdict` — Interpretive summary
- `.dds_full`, `.di_full` — Full computational objects

**When verbose=TRUE**, also returns:
- `type_s` — Sign error rate
- `type_m` — Magnitude exaggeration
- `p_sig` — P(significant)
- `dds` — Dilution Discordance Score

**dilution_scan() returns** (data.frame):
- One row per gene
- Columns: gene, bulk_beta, bulk_se, bulk_fdr, d_sc, mu_dilution, dds, di, verdict, gate_status

## Parameters

| Parameter | Default | Range | Description |
|-----------|---------|-------|-------------|
| `bulk_beta` | — | ℝ | Observed bulk effect size |
| `bulk_se` | — | (0,∞) | Standard error |
| `d_sc` | — | ℝ | Within-cluster standardized effect |
| `cluster_fraction` | — | (0,1] | f = cluster prop of tissue |
| `condition_fraction` | — | [0,1] | p = condition prev in cluster |
| `r` | 1 | (0,∞) | Relative expression (cond/tissue) |
| `alpha` | 0.05 | (0,1) | Significance level |
| `sided` | "one" | "one","two" | Test direction |
| `verbose` | FALSE | TRUE/FALSE | Include type_s, type_m, dds |

## Interpretation

### Diagnosability Index (DI)
- **< 0.30**: Underpowered. Study cannot reliably detect dilution.
- **0.30–0.80**: Uncertain power. Consider collecting more data.
- **> 0.80**: Well-powered. Study can confidently detect dilution if present.

### Dilution Discordance Score (DDS)
- **0.0–0.3**: Minimal evidence for dilution. Consistent with true null.
- **0.3–0.7**: Ambiguous. Effect could be dilution or noise.
- **0.7–1.0**: Strong to very strong evidence for dilution.

**Key caveat**: DDS is meaningless for small μ (low power). Always check DI.

## Example Results

From frozen regression test case (dpp9_m8):
```
bulk_beta = -0.084210
bulk_se = 0.062774
d_sc = 0.117962 (from single-cell)
cluster_fraction = 1332/19175 = 0.0695
condition_fraction = 0.631
alpha = 0.05, sided = "one"

Results:
  mu_dilution = 0.130536
  w = 0.0695
  di = 0.064973 (low power)
  dds = 0.454221 (weak evidence)
  type_s = 0.3685
  type_m = 15.8431
  verdict = "LOW POWER (DI=0.065): Underpowered to detect dilution"
```

Interpretation: Even though DDS suggests dilution, we lack power to conclude with confidence. Recommend collecting more cells or increasing target cluster purity.

## Regression Tests

**184 passing tests** verify:
- DDS computation (frozen: 0.454221)
- DI calculation (frozen: 0.064973)
- Type S/M error rates (frozen: 15.8431)
- r parameter across 7 values (0.5, 1, 2, 5, 10, 20, 50)
- Edge cases and error handling
- Backward compatibility with sc_zscore

Status: 0 errors, 0 warnings, 2 notes (CRAN new submission).

## Requirements

- R ≥ 4.0
- No external dependencies

## References

For complete methodology see:
- `SPEC.md` — Formal definitions and formulas
- `tests/` — Regression test specifications with frozen values
- Package vignettes for case studies

## License

[Pending: University of Washington approval]

## Authors

**Itika Arora** — Au Lab, University of Washington

---

**Version**: 0.1.0 | **Updated**: 2025-08-15 | **Repository**: https://github.com/itika05/bulknull
