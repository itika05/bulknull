test_that("DDS function computes without error on valid input", {
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
})

test_that("DDS score is always between 0 and 1", {
  scenarios <- list(
    list(
      bulk_beta = 0.0, bulk_se = 0.05, sc_zscore = 0,
      cluster_fraction = 0.05, condition_fraction = 0.5, n_cluster = 500
    ),
    list(
      bulk_beta = -0.2, bulk_se = 0.05, sc_zscore = 3.0,
      cluster_fraction = 0.05, condition_fraction = 0.5, n_cluster = 500
    ),
    list(
      bulk_beta = 0.2, bulk_se = 0.05, sc_zscore = 3.0,
      cluster_fraction = 0.05, condition_fraction = 0.5, n_cluster = 500
    )
  )

  for (scenario in scenarios) {
    result <- dilution_score(
      bulk_beta = scenario$bulk_beta,
      bulk_se = scenario$bulk_se,
      bulk_fdr = 0.5,
      sc_zscore = scenario$sc_zscore,
      cluster_fraction = scenario$cluster_fraction,
      condition_fraction = scenario$condition_fraction,
      n_cluster = scenario$n_cluster
    )

    expect_true(result$dds_score >= 0 && result$dds_score <= 1)
  }
})

test_that("DDS when cluster_fraction = 1 (entire tissue) is valid", {
  result <- dilution_score(
    bulk_beta = 0.05,
    bulk_se = 0.05,
    bulk_fdr = 0.3,
    sc_zscore = 2.0,
    cluster_fraction = 1.0,
    condition_fraction = 0.5,
    n_cluster = 1000
  )

  expect_true(result$dds_score > 0 && result$dds_score < 1)
})

test_that("DDS is low (~0.5) when single-cell z-score is 0 (no effect)", {
  result <- dilution_score(
    bulk_beta = 0.01,
    bulk_se = 0.05,
    bulk_fdr = 0.9,
    sc_zscore = 0.0,
    cluster_fraction = 0.05,
    condition_fraction = 0.5,
    n_cluster = 500
  )

  expect_true(result$dds_score < 0.6)
})
