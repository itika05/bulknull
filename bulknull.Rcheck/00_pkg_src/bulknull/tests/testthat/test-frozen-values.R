# Frozen values test: Guards against formula regressions
# Pins k=5 tissue-only values to catch any accidental changes
# Tolerance: 1e-5 for critical values; 1e-6 for derived quantities

test_that("k=5 tissue-only DPP9 values are frozen (dilution_score)", {
  result <- dilution_score(
    bulk_beta = -0.084210,
    bulk_se = 0.062774,
    bulk_fdr = 0.449745,
    sc_zscore = 2.077409,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  # Frozen k=5 tissue-only values (z-score 2.077409, tolerance 1e-5 for floating point variation)
  expect_equal(result$mu_dilution, 0.1305362, tolerance = 1e-5)
  expect_equal(result$z_bulk, -1.3414790, tolerance = 1e-5)
  expect_equal(result$dds_score, 0.4542207, tolerance = 1e-5)
})

test_that("k=5 tissue-only DPP9 DI values are frozen", {
  # Using the k=5 mu_dilution value
  mu_dilution_k5 <- 0.130536

  # One-sided test (default, recommended)
  di_one_sided <- diagnosability_index(mu_dilution = mu_dilution_k5, alpha = 0.05, sided = "one")
  expect_equal(di_one_sided$di, 0.064973, tolerance = 1e-5)

  # Two-sided test (for reference)
  di_two_sided <- diagnosability_index(mu_dilution = mu_dilution_k5, alpha = 0.05, sided = "two")
  expect_equal(di_two_sided$di, 0.033668, tolerance = 1e-5)
})

test_that("DI critical values are correct for alpha=0.05", {
  mu_test <- 0.5  # Arbitrary test value

  # One-sided: z_critical should be qnorm(1 - 0.05) = 1.644854
  di_one <- diagnosability_index(mu_dilution = mu_test, alpha = 0.05, sided = "one")
  expect_equal(di_one$critical_value, qnorm(0.95), tolerance = 1e-6)

  # Two-sided: z_critical should be qnorm(1 - 0.05/2) = 1.959964
  di_two <- diagnosability_index(mu_dilution = mu_test, alpha = 0.05, sided = "two")
  expect_equal(di_two$critical_value, qnorm(0.975), tolerance = 1e-6)
})
