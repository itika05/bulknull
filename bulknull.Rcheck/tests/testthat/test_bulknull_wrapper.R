test_that("bulknull() returns correct S3 class and structure", {
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
  expect_true("di" %in% names(result))
  expect_true("verdict" %in% names(result))
  expect_true("gate_status" %in% names(result))
})

test_that("bulknull() passes applicability gate with bulk_fdr > 0.05", {
  result <- bulknull(
    bulk_beta = -0.084210,
    bulk_se = 0.062774,
    bulk_fdr = 0.449745,
    sc_zscore = 2.077409,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_true(result$applicability$applicable)
  expect_equal(result$gate_status, "PASSED")
})

test_that("bulknull() fails applicability gate with bulk_fdr <= 0.05", {
  result <- bulknull(
    bulk_beta = -0.084210,
    bulk_se = 0.062774,
    bulk_fdr = 0.01,
    sc_zscore = 2.077409,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332,
    on_violation = "pass"
  )

  expect_false(result$applicability$applicable)
  expect_equal(result$gate_status, "FAILED")
  expect_match(result$verdict, "NOT APPLICABLE")
})

test_that("bulknull() computes correct DDS and interpretations", {
  result <- bulknull(
    bulk_beta = -0.084210,
    bulk_se = 0.062774,
    bulk_fdr = 0.449745,
    sc_zscore = 2.077409,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_true(result$dds > 0 && result$dds < 1)
  expect_match(result$dds_interpretation, "weak|Weak|Ambiguous")
})

test_that("bulknull() computes correct DI and power interpretations", {
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

  expect_true(result$di >= 0 && result$di <= 1)
  expect_match(result$di_interpretation, "Low power|Moderate power|High power")
})

test_that("bulknull() produces sensible verdicts", {
  # Realistic scenario with moderate DDS -> likely AMBIGUOUS or POSSIBLE
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

  # Low power scenario
  result_low_power <- bulknull(
    bulk_beta = -0.084210,
    bulk_se = 0.1,
    bulk_fdr = 0.2,
    sc_zscore = 0.5,
    cluster_fraction = 0.01,
    condition_fraction = 0.5,
    n_cluster = 100
  )
  expect_match(result_low_power$verdict, "UNDERPOWERED|AMBIGUOUS|DILUTION")
})

test_that("print.bulknull() executes without error", {
  result <- bulknull(
    bulk_beta = -0.084210,
    bulk_se = 0.062774,
    bulk_fdr = 0.449745,
    sc_zscore = 2.077409,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_output(print(result), "BULKNULL DIAGNOSTIC RESULT")
  expect_output(print(result), "VERDICT")
})

test_that("summary.bulknull() executes without error", {
  result <- bulknull(
    bulk_beta = -0.084210,
    bulk_se = 0.062774,
    bulk_fdr = 0.449745,
    sc_zscore = 2.077409,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_output(summary(result), "BULKNULL SUMMARY")
  expect_output(summary(result), "Status")
})

test_that("bulknull() stores full DDS and DI results internally", {
  result <- bulknull(
    bulk_beta = -0.084210,
    bulk_se = 0.062774,
    bulk_fdr = 0.449745,
    sc_zscore = 2.077409,
    cluster_fraction = 1332 / 19175,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_true(".dds_full" %in% names(result))
  expect_true(".di_full" %in% names(result))
  expect_type(result$.dds_full, "list")
  expect_type(result$.di_full, "list")
})
