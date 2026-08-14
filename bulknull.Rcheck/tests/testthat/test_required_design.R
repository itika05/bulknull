test_that("required_design() returns correct S3 class and structure", {
  result <- required_design(
    target_di = 0.3,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  expect_s3_class(result, "required_design")
  expect_type(result, "list")
  expect_true("target_di" %in% names(result))
  expect_true("required_cluster_fraction" %in% names(result))
  expect_true("required_bulk_se" %in% names(result))
  expect_true("achieved_di_at_fixed_se" %in% names(result))
  expect_true("achieved_di_at_fixed_cf" %in% names(result))
})

test_that("required_design() computes achievable targets", {
  result <- required_design(
    target_di = 0.3,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  expect_true(result$target_di == 0.3)
  expect_true(!is.na(result$required_cluster_fraction))
  expect_true(result$required_cluster_fraction > 0)
  expect_true(result$required_cluster_fraction <= 1.0)
})

test_that("required_design() handles unattainable targets with warning", {
  result <- required_design(
    target_di = 0.95,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  expect_true(is.na(result$required_cluster_fraction))
  expect_true(!is.null(result$warnings))
  expect_true(any(grepl("impossible", result$warnings, ignore.case = TRUE)))
})

test_that("required_design() supports one-sided and two-sided tests", {
  result_one <- required_design(
    target_di = 0.3,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  result_two <- required_design(
    target_di = 0.3,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "two"
  )

  expect_true(!is.na(result_one$required_cluster_fraction))
  expect_true(!is.na(result_two$required_cluster_fraction))
  # Two-sided is more conservative, so requires higher cluster fraction
  expect_true(result_two$required_cluster_fraction > result_one$required_cluster_fraction)
})

test_that("required_design() computes required bulk_se when CF fixed", {
  result <- required_design(
    target_di = 0.3,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  expect_true(result$required_bulk_se > 0)
  expect_type(result$required_bulk_se, "double")
})

test_that("print.required_design() executes without error", {
  result <- required_design(
    target_di = 0.3,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  expect_output(print(result), "Required Study Design")
  expect_output(print(result), "Target Diagnosability")
})

test_that("required_design() produces consistent interpretations", {
  result <- required_design(
    target_di = 0.5,
    sc_zscore = 2.077409,
    condition_fraction = 0.631,
    n_cluster = 1332,
    cluster_fraction = 1332 / 19175,
    bulk_se = 0.062774,
    alpha = 0.05,
    sided = "one"
  )

  expect_match(result$scenario_fixed_cf, "cluster_fraction|bulk_se")
  expect_match(result$scenario_fixed_se, "cluster_fraction|bulk_se")
})
