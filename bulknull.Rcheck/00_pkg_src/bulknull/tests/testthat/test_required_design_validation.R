# Regression tests validating required_design() against independent Python
# implementation. Input parameters:
# sc_zscore=2.077409, condition_fraction=0.631, n_cluster=1332,
# cluster_fraction=1332/19175, bulk_se=0.062774, alpha=0.05, one-sided

test_that("required_design intermediate values match Python reference", {
  sc_zscore <- 2.077409
  condition_fraction <- 0.631
  n_cluster <- 1332
  bulk_se <- 0.062774

  # Compute intermediates
  n_eff_sc <- condition_fraction * (1 - condition_fraction) * n_cluster
  d_sc <- sc_zscore / sqrt(n_eff_sc)

  expect_equal(n_eff_sc, 310.141548, tolerance = 1e-5)
  expect_equal(d_sc, 0.117962, tolerance = 1e-5)
})

test_that("required_design(target_di=0.5) matches Python", {
  result <- required_design(
    target_di = 0.5,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  expect_equal(result$required_cluster_fraction, 0.875316, tolerance = 1e-5)
  expect_equal(result$required_bulk_se, 0.0049818, tolerance = 1e-5)
})

test_that("required_design(target_di=0.8) unattainable, matches Python", {
  result <- required_design(
    target_di = 0.8,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  expect_true(is.na(result$required_cluster_fraction))
  expect_equal(result$required_bulk_se, 0.0032955417, tolerance = 1e-5)
  expect_equal(result$required_cluster_fraction_raw, 1.3231889, tolerance = 1e-6)
})

test_that("required_design(target_di=0.3) matches Python", {
  result <- required_design(
    target_di = 0.3,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  expect_equal(result$required_cluster_fraction, 0.596254, tolerance = 1e-5)
})

test_that("max DI at cluster_fraction=1.0 matches Python", {
  sc_zscore <- 2.077409
  condition_fraction <- 0.631
  n_cluster <- 1332
  bulk_se <- 0.062774
  alpha <- 0.05

  n_eff_sc <- condition_fraction * (1 - condition_fraction) * n_cluster
  d_sc <- sc_zscore / sqrt(n_eff_sc)

  max_cf <- 1.0
  d_bulk_max <- max_cf * d_sc
  mu_dilution_max <- d_bulk_max / bulk_se
  di_max <- stats::pnorm(mu_dilution_max - stats::qnorm(1 - alpha))

  expect_equal(di_max, 0.592624, tolerance = 1e-5)
})
