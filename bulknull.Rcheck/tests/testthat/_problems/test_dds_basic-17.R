# Extracted from test_dds_basic.R:17

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "bulknull", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
result <- dilution_score(
    bulk_beta = -0.03,
    bulk_se = 0.06,
    bulk_fdr = 0.82,
    sc_zscore = 2.08,
    cluster_fraction = 0.07,
    condition_fraction = 0.63,
    n_cluster = 1332
  )
expect_type(result, "list")
expect_named(
    result,
    c("dds_score", "z_bulk", "mu_dilution", "w", "d_sc", "log_lik_h1", "log_lik_h0", "summary",
      "applicable", "applicability_note", "call")
  )
