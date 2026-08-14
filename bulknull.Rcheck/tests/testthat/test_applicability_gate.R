test_that("Applicability gate correctly identifies bulk-NULL cases (applicable)", {
  result <- dilution_score(
    bulk_beta = -0.03,
    bulk_se = 0.06,
    bulk_fdr = 0.15,
    sc_zscore = 2.08,
    cluster_fraction = 0.07,
    condition_fraction = 0.63,
    n_cluster = 1332,
    on_violation = "pass"
  )

  expect_true(result$applicable)
  expect_match(result$applicability_note, "Applicable")
})

test_that("Applicability gate warns on bulk-SIGNIFICANT cases", {
  expect_warning({
    result <- dilution_score(
      bulk_beta = 3.0157,
      bulk_se = 0.6536,
      bulk_fdr = 0.0003,
      sc_zscore = 2.8455,
      cluster_fraction = 0.1492,
      condition_fraction = 0.922,
      n_cluster = 2860,
      on_violation = "warn"
    )
  }, "Not Applicable")

  # After warning, check that result is still returned
  result <- suppressWarnings({
    dilution_score(
      bulk_beta = 3.0157,
      bulk_se = 0.6536,
      bulk_fdr = 0.0003,
      sc_zscore = 2.8455,
      cluster_fraction = 0.1492,
      condition_fraction = 0.922,
      n_cluster = 2860,
      on_violation = "warn"
    )
  })
  expect_false(result$applicable)
  expect_match(result$applicability_note, "Not Applicable")
})

test_that("Applicability gate stops with error in strict mode on bulk-SIGNIFICANT", {
  expect_error({
    dilution_score(
      bulk_beta = 3.0157,
      bulk_se = 0.6536,
      bulk_fdr = 0.0003,
      sc_zscore = 2.8455,
      cluster_fraction = 0.1492,
      condition_fraction = 0.922,
      n_cluster = 2860,
      on_violation = "error"
    )
  }, "Not Applicable")
})

test_that("Applicability gate returns results silently in pass mode", {
  result <- dilution_score(
    bulk_beta = 3.0157,
    bulk_se = 0.6536,
    bulk_fdr = 0.0003,
    sc_zscore = 2.8455,
    cluster_fraction = 0.1492,
    condition_fraction = 0.922,
    n_cluster = 2860,
    on_violation = "pass"
  )

  expect_false(result$applicable)
  expect_true(is.numeric(result$dds_score))
})

test_that("Applicability gate respects custom FDR threshold", {
  result <- dilution_score(
    bulk_beta = -0.03,
    bulk_se = 0.06,
    bulk_fdr = 0.05,
    sc_zscore = 2.08,
    cluster_fraction = 0.07,
    condition_fraction = 0.63,
    n_cluster = 1332,
    null_fdr_threshold = 0.01,
    on_violation = "pass"
  )

  expect_true(result$applicable)
})

test_that("check_dds_applicability validates input FDR ranges", {
  expect_error(
    check_dds_applicability(bulk_fdr = -0.1),
    "bulk_fdr must be numeric"
  )

  expect_error(
    check_dds_applicability(bulk_fdr = 1.5),
    "bulk_fdr must be numeric"
  )

  expect_error(
    check_dds_applicability(bulk_fdr = 0.5, null_fdr_threshold = -0.05),
    "null_fdr_threshold must be numeric"
  )
})
