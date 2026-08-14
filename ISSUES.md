# Deferred Issues and Future Work

## License and Copyright Holder (PENDING)

**Status:** Pending Alfaisal University Technology Transfer Office clearance.

**Note:** The package is currently licensed under GPL-3 provisionally. The copyright holder string in DESCRIPTION is currently "Alfaisal University" (pending confirmation). The final license and copyright attribution depend on the University's IP Policy review and Technology Transfer Office approval. Do not create a public repository or DOI until this is confirmed.

---

## Within-M8 IPF-vs-Control Effect Estimation

**Status:** Deferred; not required for current claims.

**Issue:** The worked example (DPP9 in M8, k=5 tissue-only) uses a hypothesised within-M8 IPF-vs-control effect size d_sc = 0.117962, derived from the observed module enrichment z-score (M8 vs. other 11 myeloid states) as a computational anchor. This is not a measured IPF-vs-control effect.

**What is needed:** A donor-level within-M8 IPF-vs-control pseudobulk DPP9 estimate, computed as follows:
1. For each donor, compute mean DPP9 expression in M8 cells only
2. Perform donor-level limma differential expression (case vs. control donors)
3. Extract z-score from the moderated t-statistic
4. Use this z-score instead of (or in addition to) the module enrichment score

**Why:** This would allow the worked example to use a measured rather than hypothesised effect, strengthening the claims about study power and dilution diagnosis.

**Requirements:**
- Access to donor ID metadata in single-cell Seurat/SCE object
- Ability to subset by cluster and compute pseudobulk per donor
- Coordinated limma analysis at donor level

**Implementation:** Candidate function `sc_zscore_from_seurat_donor_level()` or equivalent, with mandatory `donor_col` parameter and `n_eff` computed from donor count rather than cell count.

---

## Removed / Not Implemented

### sc_zscore_from_seurat() and sc_zscore_from_sce()

**Status:** DELETED (commit: `resolution-4-implementation`)

**Reason:** These functions computed cell-level Welch's t-tests (M8 vs. other myeloid states), which do not measure IPF-vs-control effects. They were removed to prevent misinterpretation of between-state enrichment scores as within-state case-vs-control effects. Per RESOLUTION 4, sc_zscore is now interpreted only as a scale anchor for hypothesised effects.

**Future work:** If within-state case-vs-control z-scores are to be computed, implement them as donor-level pseudobulk statistics (see above), not cell-level statistics.
