# Extracted from test_bulknull_wrapper.R:99

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "bulknull", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
result <- bulknull(
    bulk_beta = -0.084210,
    bulk_se = 0.062774,
    bulk_fdr = 0.2,
    sc_zscore = 2.077409,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332
  )
expect_match(result$verdict, "DILUTION|AMBIGUOUS|UNDERPOWERED")
