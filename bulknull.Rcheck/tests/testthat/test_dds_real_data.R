test_that("Real DPP9/M8 case produces sensible DDS (0 < DDS < 1)", {
  result <- dilution_score(
    bulk_beta = -0.0297,
    bulk_se = 0.0595,
    bulk_fdr = 0.821,
    sc_zscore = 2.077,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_true(result$dds_score > 0)
  expect_true(result$dds_score < 1)
  expect_true(!is.na(result$dds_score))
  expect_true(!is.nan(result$dds_score))
})

test_that("DDS returns numeric results for various inputs", {
  result1 <- dilution_score(
    bulk_beta = 0.05,
    bulk_se = 0.1,
    bulk_fdr = 0.2,
    sc_zscore = 1.5,
    cluster_fraction = 0.15,
    condition_fraction = 0.6,
    n_cluster = 500
  )

  expect_true(is.numeric(result1$dds_score))
  expect_true(is.numeric(result1$z_bulk))
  expect_true(is.numeric(result1$mu_dilution))
})

test_that("DDS summary field contains interpretable text", {
  result <- dilution_score(
    bulk_beta = -0.03,
    bulk_se = 0.06,
    bulk_fdr = 0.82,
    sc_zscore = 2.08,
    cluster_fraction = 0.07,
    condition_fraction = 0.63,
    n_cluster = 1332
  )

  expect_true(is.character(result$summary))
  expect_true(nchar(result$summary) > 0)
  expect_true(grepl("evidence", result$summary, ignore.case = TRUE))
})
