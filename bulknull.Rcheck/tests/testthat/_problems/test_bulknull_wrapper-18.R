# Extracted from test_bulknull_wrapper.R:18

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "bulknull", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
result <- bulknull(
    bulk_beta = -0.084210,
    bulk_se = 0.062774,
    bulk_fdr = 0.449745,
    sc_zscore = 2.077409,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332,
    alpha = 0.05,
    sided = "one"
  )
expect_s3_class(result, "bulknull")
expect_type(result, "list")
expect_true("inputs" %in% names(result))
expect_true("applicability" %in% names(result))
expect_true("dds" %in% names(result))
