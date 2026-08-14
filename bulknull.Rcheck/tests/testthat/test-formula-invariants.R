# Documentation drift test: Guards against formula regressions in roxygen
# Extracts documented values and asserts they match current implementation

test_that("diagnosability_index() critical value matches documentation", {
  # From roxygen @details line 23: "Critical value: z = qnorm(1 - alpha) = 1.644854 for alpha=0.05"
  # Extract from roxygen docstring (R/diagnosability_index.R lines 22-26)
  documented_critical_one_sided <- 1.644854
  computed_critical_one_sided <- qnorm(1 - 0.05)

  expect_equal(
    computed_critical_one_sided,
    documented_critical_one_sided,
    tolerance = 1e-6,
    info = sprintf(
      "diagnosability_index() critical value mismatch: documented=%.6f, computed=%.6f (R/diagnosability_index.R line 23)",
      documented_critical_one_sided,
      computed_critical_one_sided
    )
  )

  # Verify the function returns this critical value
  di_result <- diagnosability_index(mu_dilution = 0.5, alpha = 0.05, sided = "one")
  expect_equal(
    di_result$critical_value,
    documented_critical_one_sided,
    tolerance = 1e-6,
    info = sprintf(
      "diagnosability_index(sided='one') returns wrong critical_value: %.6f (expected %.6f)",
      di_result$critical_value,
      documented_critical_one_sided
    )
  )
})

test_that("dilution_score() formula uses LINEAR dilution, not sqrt", {
  # From roxygen @details line 62 (R/dilution_score.R): "Expected bulk effect size:
  # d_bulk_expected = cluster_fraction * d_sc" (LINEAR)
  # Pre-bugfix used sqrt(cluster_fraction), which would fail this test.
  # Guard test via functional verification: use documented formula step-by-step
  # Example: k=5 tissue-only case
  sc_z <- 2.077
  cluster_frac <- 1332 / 19175
  condition_frac <- 0.631
  n_clust <- 1332
  bulk_se <- 0.062774

  # Correct formula (from roxygen)
  n_eff_sc <- condition_frac * (1 - condition_frac) * n_clust
  d_sc <- sc_z / sqrt(n_eff_sc)
  d_bulk_expected_linear <- cluster_frac * d_sc  # LINEAR
  mu_dilution_linear <- d_bulk_expected_linear / bulk_se

  # Buggy formula (for comparison, should NOT match)
  d_bulk_expected_sqrt <- sqrt(cluster_frac) * d_sc  # WRONG
  mu_dilution_sqrt <- d_bulk_expected_sqrt / bulk_se

  # The function should match the LINEAR formula
  result <- dilution_score(
    bulk_beta = -0.084210,
    bulk_se = bulk_se,
    bulk_fdr = 0.449745,
    sc_zscore = sc_z,
    cluster_fraction = cluster_frac,
    condition_fraction = condition_frac,
    n_cluster = n_clust
  )

  expect_equal(
    result$mu_dilution,
    mu_dilution_linear,
    tolerance = 1e-5,
    info = sprintf(
      "dilution_score() mu_dilution does not match LINEAR formula. Computed=%.6f, Linear=%.6f, Sqrt=%.6f (R/dilution_score.R line ~195)",
      result$mu_dilution,
      mu_dilution_linear,
      mu_dilution_sqrt
    )
  )

  expect_false(
    abs(result$mu_dilution - mu_dilution_sqrt) < 1e-5,
    info = sprintf(
      "dilution_score() matches pre-bugfix SQRT formula (%.6f), not LINEAR formula (%.6f). Regression detected.",
      mu_dilution_sqrt,
      mu_dilution_linear
    )
  )
})
