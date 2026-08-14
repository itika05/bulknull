# bulknull Statistical Specification

This document defines every statistical quantity used in the bulknull package. It records the contrasts measured at different levels and provides exact formulas, R expressions, and source provenance.

---

## RESOLUTION OF CONTRAST AND COMPOSITIONAL QUESTIONS (RESOLUTION 4)

**Background:** bulk_beta measures the IPF **sample-level** contrast (IPF case vs. control samples, n=10 vs. 7). A state-enriched effect (M8 enriched for DPP9 vs. other myeloid states) reaches bulk measurement through two potential routes:

1. **State-intrinsic route (within-M8):** A true within-M8 IPF-vs-control effect on DPP9, contributing approximately f_cluster × d_within to bulk z-score.
2. **Compositional route:** A change in M8 abundance between IPF and control, combined with M8's transcriptional distinctiveness (module enrichment score), contributing to bulk effect.

**Compositional route is closed:** The companion manuscript reports M8 abundance is **lower in IPF than control** (median 5.45% in IPF vs. 9.88% in control, bootstrap mean difference -9.2 percentage points). If M8 abundance is lower in IPF, any M8-driven effect would predict a bulk effect **opposite in direction** to the hypothesis (negative if M8 is depleted). This rules out the compositional explanation.

**State-intrinsic route is unmeasured:** A donor-level within-M8 IPF-vs-control DPP9 estimate does not exist in the project data.

**RESOLUTION 4 ADOPTED:** sc_zscore is henceforth interpreted as a **hypothesised within-state standardised effect size**, not a measured IPF-vs-control effect. Specifically:

- d_sc (derived from sc_zscore) is treated as a user-supplied order-of-magnitude anchor for a hypothesised within-M8 effect.
- In the worked example, d_sc = 0.117962 is derived from the observed module enrichment z-score (M8 vs. other myeloid) as a computational anchor point, NOT a measured case-vs-control effect.
- The framework asks: "IF there were a within-M8 effect of size d_sc, would the study design have power to detect it in bulk (DI), and how likely would the bulk-level null reflect dilution (DDS)?"
- sc_zscore_from_seurat() and sc_zscore_from_sce() are removed or replaced with donor-level pseudobulk computation, to prevent misinterpretation of cell-level scores as within-state effects.

**Interpretation:** required_design() becomes the primary entry point for study design questions, since it explicitly takes a hypothesised effect as an input parameter rather than claiming to have measured one.

---

## CORE QUANTITIES (DPP9 case study, k=5 tissue-only)

### Sample-Level Quantities

#### bulk_beta — Log Fold Change (bulk level)

**Symbol:** β_bulk

**Meaning:** Log2 fold change in DPP9 expression comparing IPF case samples to control samples, at the whole-sample (bulk) level.

**Contrast:** IPF samples vs. control samples

**Unit of analysis:** Sample

**Three variants exist in project data:**

| Variant | beta | SE | FDR | Cohort composition | Status |
|---------|------|----|----|-------------------|--------|
| **k=8 pooled** | -0.0297 | 0.0595 | 0.821 | 8 IPF + control series (includes GSE10667) | Superseded (scaling artifacts in GSE10667) |
| **k=5 tissue-only** | -0.084210 | 0.062774 | 0.449745 | 5 IPF series (GSE119600, GSE132347, GSE163033, GSE154039, GSE199926; GSE10667 excluded) | **EXAMPLE ONLY** (provenance unresolved) |
| **k=7 lung-only** | +0.0185 | unknown | unknown | 7 lung tissue series (tissue-restricted, unknown sample composition) | Not used in package |

**Status:** The k=5 tissue-only estimate is used in the worked example (dpp9_m8 data object and calculations below) for illustrative purposes only. The authoritative bulk effect estimate for this gene and cohort is NOT determined by this package. Three estimates from the project suggest different conclusions (FDR 0.821, 0.450, and unknown, depending on cohort). Which is canonical should be decided by the companion manuscript or study protocol, not by bulknull.

**Source (k=5 tissue-only, worked example):**
- File: `/Users/dritikaarora/dev/bulknull/data-raw/01_generate_dpp9_m8.R`
- Value extracted from: config/dilution_parameters.csv, row containing -0.084210
- Original derivation: limma differential expression

**Formula:** β_bulk = log2(mean_case_expression / mean_control_expression)

**R expression:**
```r
bulk_beta <- -0.084210  # k=5 tissue-only
```

**Numeric value (k=5):** -0.084210 (6 dp: -0.084210)

**Interpretation (k=5):** DPP9 is slightly lower in IPF cases than controls (fold change 2^(-0.084210) ≈ 0.944, i.e., 5.6% reduction), but this effect is not statistically significant (bulk_fdr = 0.449745).

---

#### bulk_se — Standard Error of bulk_beta

**Symbol:** SE_bulk

**Meaning:** Standard error of the log fold change estimate at the bulk level, derived from the limma linear model.

**Contrast:** Same as bulk_beta (IPF vs. control)

**Unit of analysis:** Sample

**Source:**
- File: `/Users/dritikaarora/dev/bulknull/data-raw/01_generate_dpp9_m8.R`
- Value extracted from: config/dilution_parameters.csv, row containing 0.062774
- Original derivation: limma differential expression

**Formula:** SE_bulk = sqrt(variance(β_bulk)) from limma model fit

**R expression:**
```r
bulk_se <- 0.062774
```

**Numeric value:** 0.062774 (6 dp: 0.062774)

---

#### bulk_fdr — False Discovery Rate (bulk level)

**Symbol:** FDR_bulk

**Meaning:** Adjusted p-value (FDR) from the differential expression test at the bulk level. Tests the null hypothesis that DPP9 has zero log fold change in IPF vs. control.

**Contrast:** IPF vs. control

**Unit of analysis:** Sample

**Source:**
- File: `/Users/dritikaarora/dev/bulknull/data-raw/01_generate_dpp9_m8.R`
- Value extracted from: config/dilution_parameters.csv, row containing 0.449745
- Original derivation: limma differential expression, Benjamini-Hochberg correction

**Formula:** FDR_bulk = Benjamini-Hochberg adjusted p-value from t-test or limma moderated t-statistic

**R expression:**
```r
bulk_fdr <- 0.449745
```

**Numeric value:** 0.449745 (6 dp: 0.449745)

**Interpretation:** FDR > 0.05, so bulk effect is null and framework applies.

---

### Cell-Level Quantities

#### sc_zscore — Single-Cell Z-Score (cell type module enrichment)

**Symbol:** z_sc

**Meaning:** A standardized score measuring the transcriptional distinctiveness of the M8 myeloid cluster relative to the other eleven myeloid cell states, computed from single-cell expression data. Used as a scale anchor for the hypothesised within-state effect size d_sc.

**Contrast:** M8 myeloid cluster (n=1332 cells) vs. OTHER eleven myeloid states (n=18,000+ cells, union of all non-M8 myeloid cells in k=5 tissue-only cohort)

**Unit of analysis:** Single cell (cell-type-level comparison, not sample-level)

**Source:**
- File: `/Users/dritikaarora/dev/bulknull/data-raw/01_generate_dpp9_m8.R`
- Value extracted from: config/dilution_parameters.csv, row containing 2.077409
- Original derivation: Seurat FindMarkers or equivalent (Wilcoxon rank-sum test, z-score normalization)

**Formula:**
```
z_sc = (mean_M8 - mean_other_myeloid) / pooled_se
```
where pooled_se is derived from the variance structure of the single-cell test (exact formula depends on test used; Seurat defaults to log-normalized Wilcoxon).

**R expression:**
```r
sc_zscore <- 2.077409  # M8 vs. other myeloid states (ALL samples pooled)
```

**Numeric value:** 2.077409 (6 dp: 2.077409)

**CRITICAL WARNING:** This z-score compares **M8 vs. OTHER MYELOID STATES** across all samples combined. It is NOT an IPF-vs-control comparison within M8. It is used ONLY to set the magnitude scale for d_sc in the hypothetical: "If there were a within-M8 IPF-vs-control effect with the same standardised magnitude as the M8 enrichment vs. other clusters, what power would we have?"

Do NOT interpret sc_zscore as measuring an IPF-specific effect. The framework is agnostic to the source of d_sc; it accepts it as a user input.

---

#### cluster_fraction — Cluster Size Relative to Total

**Symbol:** f_cluster

**Meaning:** The fraction of cells in the target cluster (M8) relative to all cells in the k=5 tissue-only cohort.

**Contrast:** Not a contrast; a population descriptor. M8 cells / total cells.

**Unit of analysis:** Single cell (aggregate proportion)

**Source:**
- File: `/Users/dritikaarora/dev/bulknull/data-raw/01_generate_dpp9_m8.R`
- Computed from: n_cluster = 1332, n_total_k5_tissue_only = 19175
- Value: 1332 / 19175

**Formula:** f_cluster = n_cluster / n_total

**R expression:**
```r
cluster_fraction <- 1332 / 19175
# = 0.069465451...
```

**Numeric value:** 0.069465 (6 dp: 0.069465)

---

#### condition_fraction — Case Fraction within Cohort

**Symbol:** f_condition

**Meaning:** The fraction of cells in the case (IPF) condition, across all cells in the k=5 tissue-only cohort.

**Contrast:** IPF case cells / total cells (across all clusters and samples)

**Unit of analysis:** Single cell (aggregate proportion)

**Source:**
- File: `/Users/dritikaarora/dev/bulknull/data-raw/01_generate_dpp9_m8.R`
- Value extracted from: config/dilution_parameters.csv, row containing 0.631
- Original derivation: Cell count from k=5 tissue-only metadata: (10 IPF samples * avg cells per IPF sample) / total cells across 10 IPF + 7 control samples
- Justification: 5 IPF series (GSE119600, GSE132347, GSE163033, GSE154039, GSE199926) have case/control sample metadata

**Formula:** f_condition = (sum of cells in all case samples) / (total cells across all samples in cohort)

**R expression:**
```r
condition_fraction <- 0.631
```

**Numeric value:** 0.631000 (6 dp: 0.631000)

---

#### n_cluster — Cell Count in Target Cluster

**Symbol:** n_M8

**Meaning:** The total number of cells assigned to the M8 myeloid cluster in the k=5 tissue-only cohort.

**Contrast:** Not a contrast; a count.

**Unit of analysis:** Single cell

**Source:**
- File: `/Users/dritikaarora/dev/bulknull/data-raw/01_generate_dpp9_m8.R`
- Value extracted from: config/dilution_parameters.csv, row containing 1332
- Original derivation: Cluster assignments from k=5 tissue-only Seurat object

**Formula:** n_cluster = sum(cluster_assignment == "M8")

**R expression:**
```r
n_cluster <- 1332
```

**Numeric value:** 1332 (exact integer)

---

## DERIVED QUANTITIES

### n_eff — Effective Sample Size (single-cell)

**Symbol:** n_eff_sc

**Meaning:** The effective sample size for the single-cell test, accounting for the case-control balance within the cohort. Used to scale the single-cell effect size to a per-subject metric.

**Formula:**
```
n_eff_sc = f_condition * (1 - f_condition) * n_cluster
```

where:
- f_condition = 0.631
- (1 - f_condition) = 0.369
- n_cluster = 1332

**R expression:**
```r
n_eff_sc <- condition_fraction * (1 - condition_fraction) * n_cluster
         <- 0.631 * 0.369 * 1332
         <- 310.141548
```

**Numeric value:** 310.141548 (6 dp: 310.141548)

**Interpretation:** The case-control imbalance (0.631 vs 0.369) and cluster size together give an effective sample size of ~310 for the single-cell signal.

---

### d_sc — Hypothesised Within-State Standardised Effect Size

**Symbol:** d_sc

**Meaning:** A hypothesised standardised effect size for DPP9 within the M8 myeloid cluster, in IPF vs. control, expressed per unit effective sample size. This is NOT a measured IPF-vs-control effect. Instead, it is supplied as an order-of-magnitude anchor for counterfactual study design questions: "IF there were a within-M8 effect of size d_sc, what power would the study have?"

**Status:** Derived from the observed module enrichment z-score (M8 vs. other myeloid states) as a computational anchor, NOT measured as an IPF-vs-control effect within M8.

**Formula (computational):**
```
d_sc = z_sc / sqrt(n_eff_sc)
     = 2.077409 / sqrt(310.141548)
```

where z_sc is the between-state enrichment score (M8 vs. other 11 myeloid states) and n_eff_sc is the effective sample size accounting for case-control imbalance.

**R expression:**
```r
d_sc <- sc_zscore / sqrt(n_eff_sc)
     <- 2.077409 / sqrt(310.141548)
     <- 0.117962
```

**Numeric value:** 0.117962 (6 dp: 0.117962)

**Important caveat:** This value is an ANCHOR derived from observed M8 enrichment (between-state comparison), not a measured within-M8 IPF-vs-control effect. The framework uses d_sc to ask "what if there were a within-state effect of this magnitude?" — not "there is a within-state effect of this magnitude."

---

### mu_dilution — Expected Bulk Effect under Dilution Hypothesis

**Symbol:** μ_dilution

**Meaning:** The expected log fold change (in bulk measurement) if the bulk-level null (β_bulk ≈ 0) is explained by dilution: the true effect exists in M8 but is diluted across the whole sample because M8 is only 6.95% of cells.

**Formula:**
```
μ_dilution = (f_cluster * d_sc) / bulk_se
           = (0.069465 * 0.117962) / 0.062774
```

**R expression:**
```r
mu_dilution <- (cluster_fraction * d_sc) / bulk_se
            <- (0.069465 * 0.117962) / 0.062774
            <- 0.130536
```

**Numeric value:** 0.130536 (6 dp: 0.130536)

**Interpretation:** If the true effect is d_sc per unit in M8, and M8 is 6.95% of cells, the bulk measurement would see z_bulk ≈ 0.1305 under the dilution model.

---

### z_bulk — Observed Bulk Z-Score

**Symbol:** z_bulk

**Meaning:** The standard normal z-score corresponding to the observed bulk log fold change, i.e., how many standard errors away from zero the bulk estimate is.

**Formula:**
```
z_bulk = beta_bulk / se_bulk
       = -0.084210 / 0.062774
```

**R expression:**
```r
z_bulk <- bulk_beta / bulk_se
       <- -0.084210 / 0.062774
       <- -1.341479
```

**Numeric value:** -1.341479 (6 dp: -1.341479)

**Interpretation:** The observed bulk effect is 1.34 standard errors below zero, consistent with a true null.

---

### DDS — Dilution Discordance Score

**Symbol:** DDS

**Meaning:** The Bayesian posterior probability that the observed bulk z-score came from the dilution hypothesis (true effect in M8, diluted in bulk) rather than the genuine null hypothesis (true zero effect everywhere).

**Contrasts:** Compares two models:
- H1 (dilution): z_bulk ~ N(μ_dilution, 1)
- H0 (null): z_bulk ~ N(0, 1)

**Formula:**
```
DDS = φ(z_bulk - μ_dilution) / [φ(z_bulk - μ_dilution) + φ(z_bulk)]
```

where φ(x) is the standard normal PDF: φ(x) = (1/√(2π)) * exp(-x²/2)

**R expression:**
```r
mu_dilution <- 0.130536
z_bulk <- -1.341479

numerator <- stats::dnorm(z_bulk - mu_dilution)
          <- stats::dnorm(-1.341479 - 0.130536)
          <- stats::dnorm(-1.472015)
          <- 0.1381...

denominator <- numerator + stats::dnorm(z_bulk)
            <- 0.1381... + stats::dnorm(-1.341479)
            <- 0.1381... + 0.1758...
            <- 0.3139...

DDS <- numerator / denominator
    <- 0.4542207
```

**Numeric value:** 0.454221 (6 dp: 0.454221)

**Interpretation:** ~45.4% posterior probability of dilution vs. ~54.6% posterior probability of genuine null.

---

### DI — Diagnosability Index

**Symbol:** DI

**Meaning:** Statistical power to detect dilution if it truly occurred at the observed single-cell effect size (d_sc), given the study design. The probability that a one-sided Welch's t-test would reject the null at significance level α, under the alternative that the true bulk effect is μ_dilution.

**Formula (one-sided, α = 0.05):**
```
z_critical = Φ^{-1}(1 - α) = qnorm(1 - 0.05) = 1.644854

DI = Φ(μ_dilution - z_critical)
   = Φ(0.130536 - 1.644854)
   = Φ(-1.514318)
```

where Φ(x) is the standard normal CDF.

**R expression:**
```r
z_critical <- stats::qnorm(1 - 0.05)
           <- 1.644854

di_one <- stats::pnorm(mu_dilution - z_critical)
       <- stats::pnorm(0.130536 - 1.644854)
       <- stats::pnorm(-1.514318)
       <- 0.064973
```

**Numeric value (one-sided):** 0.064973 (6 dp: 0.064973)

**Formula (two-sided, α = 0.05):**
```
z_critical_two = Φ^{-1}(1 - α/2) = qnorm(0.975) = 1.959964

DI_two = Φ(μ_dilution - z_critical_two)
       = Φ(0.130536 - 1.959964)
       = Φ(-1.829428)
```

**R expression:**
```r
z_critical_two <- stats::qnorm(1 - 0.05/2)
               <- 1.959964

di_two <- stats::pnorm(mu_dilution - z_critical_two)
       <- stats::pnorm(0.130536 - 1.959964)
       <- stats::pnorm(-1.829428)
       <- 0.033668
```

**Numeric value (two-sided):** 0.033668 (6 dp: 0.033668)

**Interpretation:** Only 6.5% power (one-sided) to detect dilution at this effect size. The study is severely underpowered.

---

## FUNCTION ARGUMENTS AND RETURN VALUES

### bulknull()

**Purpose:** One-call wrapper that computes applicability gate, DDS, and DI, returning an S3 object with interpretations.

**Arguments:**

| Argument | Type | Symbol | Definition | Example Value |
|----------|------|--------|-----------|----------------|
| bulk_beta | numeric | β_bulk | Log fold change from bulk DE test | -0.084210 |
| bulk_se | numeric | SE_bulk | Standard error of bulk_beta | 0.062774 |
| bulk_fdr | numeric | FDR_bulk | FDR from bulk DE test; applicability gate checks if > null_fdr_threshold | 0.449745 |
| sc_zscore | numeric | z_sc | Standardized effect size from single-cell test (M8 vs. other myeloid) | 2.077409 |
| cluster_fraction | numeric | f_cluster | M8 cells / total k=5 cells | 0.069465 |
| condition_fraction | numeric | f_condition | IPF case cells / total k=5 cells | 0.631 |
| n_cluster | numeric | n_M8 | Total cells in M8 | 1332 |
| null_fdr_threshold | numeric | α_gate | FDR threshold for applicability gate; default 0.05 | 0.05 |
| on_violation | character | — | Behavior if bulk_fdr ≤ null_fdr_threshold: "error" (stop), "warn" (warn + return NA), "pass" (return result marked FAILED) | "error" |
| alpha | numeric | α_di | Significance level for DI power calculation; default 0.05 | 0.05 |
| sided | character | — | One-sided or two-sided test for DI; default "one" | "one" |

**Returns:** S3 object of class "bulknull" containing:
- inputs (list): all input parameters
- applicability (logical): whether framework applies (bulk_fdr > threshold)
- dds (numeric): Dilution Discordance Score, value 0–1
- mu_dilution (numeric): expected bulk z-score under dilution
- z_bulk (numeric): observed bulk z-score
- dds_interpretation (character): categorical band ("Strong", "Weak-to-moderate", "Weak", "Minimal")
- di (numeric): Diagnosability Index, value 0–1
- di_interpretation (character): power band ("Low", "Moderate", "High")
- verdict (character): summary decision string
- gate_status (character): "PASSED" or "FAILED"

---

### required_design()

**Purpose:** Inverts the DI formula to find study parameters needed to achieve a target power level.

**Arguments:**

| Argument | Type | Symbol | Definition | Example Value |
|----------|------|--------|-----------|----------------|
| target_di | numeric | DI_target | Desired diagnosability index (0–1) | 0.5 |
| sc_zscore | numeric | z_sc | Single-cell z-score (same as bulknull) | 2.077409 |
| condition_fraction | numeric | f_condition | Case fraction (same as bulknull) | 0.631 |
| n_cluster | numeric | n_M8 | Cluster size (same as bulknull) | 1332 |
| cluster_fraction | numeric | f_cluster | Current cluster fraction; used in fixed-SE scenario only | 0.069465 |
| bulk_se | numeric | SE_bulk | Current bulk SE; used as fixed value in one scenario | 0.062774 |
| alpha | numeric | α_di | Significance level; default 0.05 | 0.05 |
| sided | character | — | One-sided or two-sided; default "one" | "one" |

**Scenario 1 (Solve for cluster_fraction at fixed bulk_se):**

**Formula:**
```
Solve: DI_target = Φ(μ_required - z_critical)
for:   μ_required = Φ^{-1}(DI_target) + z_critical

Then: μ_required = (f_cluster_required * d_sc) / bulk_se

Therefore: f_cluster_required = μ_required * bulk_se / d_sc
```

**Example (target_di = 0.5):**
```r
z_critical <- qnorm(1 - 0.05) = 1.644854
mu_required <- qnorm(0.5) + z_critical = 0 + 1.644854 = 1.644854
d_sc <- 0.117962
bulk_se <- 0.062774

f_cluster_required <- 1.644854 * 0.062774 / 0.117962
                   <- 0.875316
```

**Numeric value:** 0.875316 (6 dp: 0.875316) — unattainable (exceeds 1.0 at target_di > ~0.8)

**Scenario 2 (Solve for bulk_se at fixed cluster_fraction):**

**Formula:**
```
Given: f_cluster (passed in)
Compute: d_bulk_expected = f_cluster * d_sc
Then: SE_bulk_required = d_bulk_expected / μ_required
```

**Example (target_di = 0.5, f_cluster = 0.069465):**
```r
d_bulk_expected <- 0.069465 * 0.117962 = 0.008198
mu_required <- 1.644854

bulk_se_required <- 0.008198 / 1.644854 = 0.0049818
```

**Numeric value:** 0.0049818 (6 dp: 0.0049818)

**Returns:** S3 object of class "required_design" containing:
- target_di (numeric): requested power level
- required_cluster_fraction (numeric): cluster fraction needed at fixed bulk_se; NA if >1.0 (unattainable)
- required_cluster_fraction_raw (numeric): raw computed value (even if unattainable); e.g., 1.323189
- achieved_di_at_fixed_se (numeric): actual DI achieved at required_cluster_fraction
- required_bulk_se (numeric): bulk SE needed at fixed cluster_fraction
- achieved_di_at_fixed_cf (numeric): actual DI achieved at required_bulk_se
- scenario_fixed_cf (character): interpretive message for Scenario 1
- scenario_fixed_se (character): interpretive message for Scenario 2
- warnings (character vector or NULL): warnings if targets unattainable

---

### dilution_scan()

**Purpose:** Apply bulknull framework vectorized across multiple genes in a DEG table.

**Arguments:**

| Argument | Type | Definition | Example |
|----------|------|-----------|---------|
| bulk_deg | data.frame | Differential expression results with columns: gene (character), beta or log2fc (numeric effect), se or se_log2fc (numeric error), fdr or padj (numeric adjusted p-value) | data.frame(gene=c("GENE1","GENE2"), beta=c(-0.084, 0.020), se=c(0.063, 0.050), fdr=c(0.45, 0.20)) |
| sc_zscores | named numeric vector | z-scores keyed by gene name | c(GENE1=2.077, GENE2=1.5) |
| cluster_fraction | numeric | Same f_cluster for all genes | 0.069465 |
| condition_fraction | numeric | Same f_condition for all genes | 0.631 |
| n_cluster | numeric | Same n_cluster for all genes | 1332 |
| null_fdr_threshold | numeric | Applicability gate threshold; default 0.05 | 0.05 |
| alpha | numeric | DI significance level; default 0.05 | 0.05 |
| sided | character | One-sided or two-sided DI; default "one" | "one" |

**Processing:**
- For each gene in bulk_deg:
  - If gene not in sc_zscores: skip, mark as "SKIPPED"
  - If sc_zscore == 0: skip, mark as "SKIPPED"
  - If bulk_fdr ≤ null_fdr_threshold: mark gate_status = "FAILED"
  - Else: call bulknull() internally, extract results

**Returns:** data.frame with columns:
- gene (character): gene name
- bulk_beta, bulk_se, bulk_fdr (numeric): from bulk_deg
- sc_zscore (numeric): from sc_zscores or NA
- mu_dilution, dds, di (numeric): computed via bulknull()
- dds_interpretation, di_interpretation (character): categorical bands
- verdict (character): summary decision
- gate_status (character): "PASSED", "FAILED", or "SKIPPED"

---

## SUMMARY TABLE

| Quantity | Symbol | Value (6 dp) | Unit | Formula/Source |
|----------|--------|-------|------|---------|
| bulk_beta | β_bulk | -0.084210 | log2(FC) | limma, k=5 tissue-only |
| bulk_se | SE_bulk | 0.062774 | log2(FC) | limma, k=5 tissue-only |
| bulk_fdr | FDR_bulk | 0.449745 | (unitless) | Benjamini-Hochberg, limma |
| sc_zscore | z_sc | 2.077409 | (unitless) | Wilcoxon, M8 vs. other myeloid |
| cluster_fraction | f_cluster | 0.069465 | (unitless) | 1332 / 19175 |
| condition_fraction | f_condition | 0.631000 | (unitless) | case cells / total cells, k=5 |
| n_cluster | n_M8 | 1332 | (integer) | count, M8 cluster |
| n_eff | n_eff_sc | 310.141548 | (unitless) | 0.631 * 0.369 * 1332 |
| d_sc | d_sc | 0.117962 | (unitless) | 2.077409 / √310.141548 |
| mu_dilution | μ_dilution | 0.130536 | z-score | (0.069465 * 0.117962) / 0.062774 |
| z_bulk | z_bulk | -1.341479 | z-score | -0.084210 / 0.062774 |
| DDS | DDS | 0.454221 | probability | φ(z_bulk - μ) / [φ(z_bulk - μ) + φ(z_bulk)] |
| DI (one-sided) | DI_1 | 0.064973 | probability | Φ(0.130536 - 1.644854) |
| DI (two-sided) | DI_2 | 0.033668 | probability | Φ(0.130536 - 1.959964) |

---

## RESOLVED ISSUES

1. **Contrast discrepancy and compositional question (RESOLUTION 4):** bulk_beta measures IPF vs. control; sc_zscore measures M8 vs. other myeloid (between-state enrichment). Compositional route is closed (M8 lower in IPF, opposite direction). State-intrinsic route is unmeasured. d_sc is now defined as a hypothesised effect, not a measured one. required_design() is the primary entry point for study design.

## DEFERRED WORK

1. **Within-M8 IPF-vs-control effect estimation:** Computing a donor-level within-M8 IPF-vs-control DPP9 estimate would allow the worked example to use a measured rather than hypothesised effect. This requires access to donor-level metadata in the single-cell object and pseudobulk computation. Deferred; not required for current claims.

---

END OF SPEC.md
