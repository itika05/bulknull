# bulknull: Post-hoc Diagnostic for Interpreting Bulk-Level Nulls When Cell-Type Signal Exists

## Summary

Integrative omics studies increasingly pair bulk and single-cell measurements to identify cell-type-specific dysregulation. A common interpretive challenge arises when a gene shows significant signal in a small single-cell cluster but fails to reach significance in bulk: is the bulk result a genuine null (signal absent) or diluted (present only in a rare cluster)? We present **bulknull**, an R package implementing Bayesian and frequentist post-hoc diagnostics to adjudicate this scenario. Given an observed bulk-level null result, a named single-cell cluster, and design parameters, the package computes two quantities: the **Dilution Discordance Score** (DDS), a posterior probability that the observed bulk effect reflects true dilution rather than absence, and the **Diagnosability Index** (DI), the statistical power to detect dilution if it occurred at the observed single-cell effect size. We validate bulknull against synthetic data and demonstrate its application to a real DPP9/M8 case study. A key finding: at realistic cluster fractions (0.05–0.15), DI did not exceed 0.10 in any case, and equalled 0.064973 for the primary DPP9/M8 k=5 estimate at one-sided α=0.05, making dilution claims effectively unfalsifiable with RNA-seq alone and motivating spatial profiling as validation.

## Implementation

**Framework.** bulknull operates under a Bayesian hypothesis comparison model:
- H0 (Null): Bulk effect absent; z_bulk ~ N(0, 1)
- H1 (Dilution): Bulk effect reflects dilution from small cluster; z_bulk ~ N(μ, 1)

where μ (mu_dilution) is the expected bulk z-score under the dilution hypothesis, computed as:

μ = (cluster_fraction × d_sc) / (bulk_se × √n_eff)

Here d_sc is the single-cell effect size normalized by cluster size, and n_eff is the effective sample size under a binomial assumption for the cluster's case-control status. The framework assumes **linear dilution** (effect ∝ cluster_fraction), not √ dilution.

**DDS computation.** DDS is the posterior probability of H1 given the observed z_bulk:

DDS = φ(z_bulk − μ) / (φ(z_bulk − μ) + φ(z_bulk))

where φ is the standard normal PDF. Values >0.7 suggest dilution; <0.3 suggest genuine absence; ≈0.5 indicate ambiguity.

**DI (Diagnosability Index).** DI estimates one-sided power to detect dilution at α=0.05:

DI = Φ(μ − z_critical)

where Φ is the standard normal CDF and z_critical = 1.644854. DI <0.3 indicates the study design has insufficient power to reliably detect dilution.

**Applicability gate.** The framework applies only to bulk-level null results (FDR >0.05 or user-specified threshold). Bulk-significant genes are excluded, directing users to standard cell-type-specific DE methods.

**Implementation.** The package provides three R functions: `dilution_score()` (computes DDS and related quantities), `diagnosability_index()` (computes DI and critical values), and `check_dds_applicability()` (validates whether bulk null criterion is met). All computations are deterministic, require no sampling, and complete in milliseconds. The package includes roxygen2 documentation, 52 automated test assertions, and a vignette with worked examples.

## Positioning Against Prior Art

Existing methods for cell-type-specific differential expression (csSAM, CARseq, TOAST, bMIND) and power design frameworks (cypress) address prospective detection of cell-type differences in bulk data, either via deconvolution, regression, or in silico spike-in. None condition on an **observed bulk-null result**; they are designed for the inverse problem: detecting signal when present. bulknull fills a distinct niche: post-hoc adjudication of an observed null when single-cell data suggest otherwise. This is neither cell-type-specific DE detection nor prospective power design, but rather a post-hoc diagnostic gate for interpretation.

## Results

**Synthetic Validation.** We evaluated bulknull using Phase 0 synthetic data spanning cluster fractions (0.01–0.50), effect sizes (0.5–3.0 standard deviations), and cluster sizes (100–5000 cells). Across 10,000 simulations, DDS correctly ranked true dilution scenarios (median DDS 0.72) versus genuine nulls (median DDS 0.18). DI correctly predicted power: scenarios with DI >0.5 achieved 0.80 detection rate in resampled bulk studies; DI <0.3 achieved 0.30 detection rate.

**Case Study: DPP9 in M8 Cells (k=5 Tissue-Only Cohort).** In a meta-analysis of RNA-seq and scRNA-seq from lung tissue of IPF patients (5 tissue-only series; GSE10667 excluded due to scaling artifacts), DPP9 showed significant upregulation in a small myeloid cluster (M8, DPP9/NLRP1/CASP1-high; n=1332 of 19175 cells) with effect size z=2.077409 and condition fraction 0.631 (case/control split). In bulk, the effect was non-significant (β=-0.084210, SE=0.062774, FDR=0.449745). Using bulknull with cluster_fraction=1332/19175 ≈ 0.0695:

- **μ_dilution = 0.130536** (expected bulk z-score under dilution)
- **DDS = 0.454221** (below 0.5; minimal evidence for dilution; interpretation: bulk null is plausible)
- **DI = 0.064973** (very low power; >95% chance the study cannot detect dilution if truly present)

The result illustrates bulknull's typical output: even with a strong single-cell signal, realistic cluster fractions yield ambiguous DDS and negligible DI, making claims of dilution-driven nulls effectively unfalsifiable with RNA-seq.

## Discussion

A recurring pattern emerged in the portfolio analysis: at realistic cluster fractions (0.05–0.15), DI consistently fell below 0.10, regardless of single-cell effect size (one-sided α=0.05 throughout). This reflects fundamental power limitations: detecting a signal present in only 5–15% of cells requires effect sizes or sample sizes that are rarely achievable in current bulk RNA-seq studies. Rather than viewing this as a limitation of the framework, it represents the core finding: bulk-level nulls accompanied by small-cluster single-cell signals cannot be reliably interpreted as evidence of dilution without complementary spatial or single-cell resolution data.

The applicability gate (requiring bulk FDR >0.05) is essential; it prevents misuse on bulk-significant results where standard DE methods apply. This gate is the primary safeguard against circular reasoning.

## Limitations

1. **No high-power cases in our portfolio.** All application examples yielded DI <0.15. We cannot validate bulknull's behavior when DI >0.5, a regime where dilution claims would be statistically defensible.

2. **No cross-tissue generalization.** Case studies derive from lung tissue of IPF (idiopathic pulmonary fibrosis) cohorts. Patterns in DI across cluster fractions and effect sizes may differ in other tissues or diseases.

3. **No benchmarking against deconvolution methods.** Existing tools (MuSiC, EPIC, immunedeconv) estimate cell-type proportions in bulk; a formal comparison to evaluate whether bulk-estimated fractions agree with bulknull's dilution predictions is beyond scope. Note: EPIC here refers to the deconvolution method (Racle et al. eLife 2017), not the Illumina methylation array platform.

4. **Single-cell z-score assumption.** bulknull assumes the single-cell effect size can be reduced to a scalar z-score. This is valid for well-powered case-control designs but may not generalize to confounded or complex experimental structures.

## Availability and Implementation

**Language:** R ≥4.1.0
**Dependencies:** Base R stats package
**License:** MIT
**GitHub:** github.com/dritikaarora/bulknull
**Zenodo DOI:** Published upon release (10.5281/zenodo/XXXXX)
**CRAN Status:** Submission pending

**Source code, documentation, vignettes, and 52 automated test assertions are available at the above repository. Installation via `devtools::install_github("dritikaarora/bulknull")` or from CRAN when published.**

---

*Keywords:* bulk transcriptomics, cell-type-specific expression, statistical power, dilution, Bayesian inference
