# Extracted from test_dilution_scan.R:111

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "bulknull", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
bulk_deg <- data.frame(
    gene = c("DPP9"),
    beta = c(-0.084210),
    se = c(0.062774),
    fdr = c(0.449745)
  )
sc_zscores <- c(DPP9 = 2.077409)
result <- dilution_scan(
    bulk_deg = bulk_deg,
    sc_zscores = sc_zscores,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332,
    alpha = 0.05,
    sided = "one"
  )
expect_equal(result$gate_status, "PASSED")
expect_true(result$dds > 0 && result$dds < 1)
