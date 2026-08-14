test_that("DDS rejects invalid inputs with clear error messages", {
  expect_error(
    dilution_score(
      bulk_beta = 0.1,
      bulk_se = -0.05,
      bulk_fdr = 0.5,
      sc_zscore = 2.0,
      cluster_fraction = 0.05,
      condition_fraction = 0.5,
      n_cluster = 500
    ),
    "bulk_se must be positive"
  )

  expect_error(
    dilution_score(
      bulk_beta = 0.1,
      bulk_se = 0.05,
      bulk_fdr = 0.5,
      sc_zscore = 2.0,
      cluster_fraction = 1.5,
      condition_fraction = 0.5,
      n_cluster = 500
    ),
    "cluster_fraction must be in"
  )

  expect_error(
    dilution_score(
      bulk_beta = 0.1,
      bulk_se = 0.05,
      bulk_fdr = 0.5,
      sc_zscore = 2.0,
      cluster_fraction = 0.05,
      condition_fraction = 1.5,
      n_cluster = 500
    ),
    "condition_fraction must be in"
  )

  expect_error(
    dilution_score(
      bulk_beta = 0.1,
      bulk_se = 0.05,
      bulk_fdr = 0.5,
      sc_zscore = 2.0,
      cluster_fraction = 0.05,
      condition_fraction = 0.5,
      n_cluster = 0
    ),
    "n_cluster must be"
  )
})

test_that("Bulk z-score is correctly computed as beta / se", {
  beta <- 0.1234
  se <- 0.0567

  result <- dilution_score(
    bulk_beta = beta,
    bulk_se = se,
    bulk_fdr = 0.5,
    sc_zscore = 2.0,
    cluster_fraction = 0.05,
    condition_fraction = 0.5,
    n_cluster = 500
  )

  expected_z <- beta / se
  expect_equal(result$z_bulk, expected_z, tolerance = 1e-6)
})

test_that("Expected dilution z (mu_dilution) follows corrected formula", {
  sc_z <- 2.5
  cluster_frac <- 0.08
  condition_frac <- 0.6
  n_clust <- 1000
  bulk_se_test <- 0.05

  result <- dilution_score(
    bulk_beta = 0.0,
    bulk_se = bulk_se_test,
    bulk_fdr = 0.5,
    sc_zscore = sc_z,
    cluster_fraction = cluster_frac,
    condition_fraction = condition_frac,
    n_cluster = n_clust
  )

  n_eff_sc <- condition_frac * (1 - condition_frac) * n_clust
  d_sc <- sc_z / sqrt(n_eff_sc)
  d_bulk_exp <- cluster_frac * d_sc
  expected_mu <- d_bulk_exp / bulk_se_test

  expect_equal(result$mu_dilution, expected_mu, tolerance = 1e-6)
})

test_that("Results are deterministic (no randomness in DDS function)", {
  result1 <- dilution_score(
    bulk_beta = -0.03,
    bulk_se = 0.06,
    bulk_fdr = 0.82,
    sc_zscore = 2.08,
    cluster_fraction = 0.07,
    condition_fraction = 0.63,
    n_cluster = 1332
  )

  result2 <- dilution_score(
    bulk_beta = -0.03,
    bulk_se = 0.06,
    bulk_fdr = 0.82,
    sc_zscore = 2.08,
    cluster_fraction = 0.07,
    condition_fraction = 0.63,
    n_cluster = 1332
  )

  expect_equal(result1$dds_score, result2$dds_score)
  expect_equal(result1$mu_dilution, result2$mu_dilution)
})
