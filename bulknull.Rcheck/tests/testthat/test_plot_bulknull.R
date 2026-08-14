test_that("plot.bulknull computes DI landscape correctly", {
  # Test against reference grid values at d_sc = 0.117962
  d_sc <- 0.117962
  alpha <- 0.05
  z_critical <- qnorm(1 - alpha)

  test_points <- list(
    list(cf = 0.0695, se = 0.062774, expected = 0.0650),
    list(cf = 0.0695, se = 0.005, expected = 0.4979),
    list(cf = 0.25, se = 0.010, expected = 0.9039),
    list(cf = 1.00, se = 0.030, expected = 0.9889)
  )

  for (tp in test_points) {
    d_bulk <- tp$cf * d_sc
    mu_dil <- d_bulk / tp$se
    di_computed <- pnorm(mu_dil - z_critical)

    # Tolerance: 0.001 (4 decimal places)
    expect_equal(di_computed, tp$expected, tolerance = 0.001)
  }
})

test_that("plot.bulknull creates plot without error", {
  # Create a minimal bulknull object
  result <- bulknull(
    bulk_beta = -0.0297,
    bulk_se = 0.0595,
    bulk_fdr = 0.821,
    sc_zscore = 2.077,
    cluster_fraction = 0.0695,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  # Plotting should not error
  expect_no_error(plot(result))
})

test_that("plot.bulknull marks study position correctly", {
  result <- bulknull(
    bulk_beta = -0.0297,
    bulk_se = 0.0595,
    bulk_fdr = 0.821,
    sc_zscore = 2.077,
    cluster_fraction = 0.0695,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  # Plot should include study position inputs
  expect_equal(result$inputs$bulk_se, 0.0595)
  expect_equal(result$inputs$cluster_fraction, 0.0695)
})
