# Regression tests for r (relative expression) parameter
# TASK 1: Verify computed values at specified r levels
# Test parameters: f=1332/19175, d_sc=0.117962, bulk_se=0.062774, alpha=0.05 one-sided

test_that("r parameter: w and mu_dilution computed correctly at r values", {

  # Reference parameters
  d_sc <- 0.117962
  cluster_fraction <- 1332 / 19175
  bulk_se <- 0.062774
  condition_fraction <- 0.631
  n_cluster <- 1332

  # Precomputed regression values (tolerance 1e-3 per TASK 1)
  expected_values <- data.frame(
    r = c(0.5, 1.0, 2, 5, 10, 20, 50),
    w = c(0.0360, 0.0695, 0.1299, 0.2718, 0.4274, 0.5989, 0.7887),
    mu = c(0.0676, 0.1305, 0.2441, 0.5108, 0.8032, 1.1254, 1.4821),
    di = c(0.0574, 0.0650, 0.0806, 0.1284, 0.2000, 0.3017, 0.4354)
  )

  r_values <- c(0.5, 1.0, 2, 5, 10, 20, 50)

  for (i in seq_along(r_values)) {
    r <- r_values[i]
    exp <- expected_values[i, ]

    result <- dilution_score(
      bulk_beta = -0.084210,
      bulk_se = bulk_se,
      bulk_fdr = 0.449745,
      sc_zscore = 2.077409,
      cluster_fraction = cluster_fraction,
      condition_fraction = condition_fraction,
      r = r,
      n_cluster = n_cluster
    )

    expect_equal(result$w, exp[["w"]], tolerance = 1e-3,
                 label = sprintf("w at r=%g", r))
    expect_equal(result$mu_dilution, exp[["mu"]], tolerance = 1e-3,
                 label = sprintf("mu_dilution at r=%g", r))

    di_result <- diagnosability_index(
      mu_dilution = result$mu_dilution,
      alpha = 0.05,
      sided = "one"
    )

    expect_equal(di_result$di, exp[["di"]], tolerance = 1e-3,
                 label = sprintf("DI at r=%g", r))
  }
})

test_that("critical_precision() at r parameter reports unattainable at target_di=0.8", {

  # At all r values, target DI 0.8 should be unattainable
  # because mu_max (when w->1) = d_sc/bulk_se = 1.8792 < mu_required(0.8)

  d_sc <- 0.117962
  bulk_se <- 0.062774

  mu_max <- d_sc / bulk_se  # 1.8792 (when w approaches 1)

  result <- critical_precision(
    d_sc = d_sc,
    bulk_se = bulk_se,
    target_di = 0.8,
    alpha = 0.05,
    sided = "one"
  )

  # mu_required for DI 0.8: qnorm(0.8) + qnorm(0.95) = 0.8416 + 1.6449 = 2.4865
  expect_false(result$attainable,
               label = "target_di=0.8 should be unattainable (2.4865 > 1.8792)")
})
