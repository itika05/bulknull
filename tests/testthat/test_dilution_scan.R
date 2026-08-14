test_that("dilution_scan() returns data frame with correct columns", {
  bulk_deg <- data.frame(
    gene = c("GENE1", "GENE2"),
    beta = c(-0.084, 0.020),
    se = c(0.063, 0.050),
    fdr = c(0.450, 0.200)
  )
  sc_zscores <- c(GENE1 = 2.077, GENE2 = 1.5)

  result <- dilution_scan(
    bulk_deg = bulk_deg,
    sc_zscores = sc_zscores,
    cluster_fraction = 0.0695,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_s3_class(result, "data.frame")
  expect_true("gene" %in% names(result))
  expect_true("dds" %in% names(result))
  expect_true("di" %in% names(result))
  expect_true("verdict" %in% names(result))
  expect_true("gate_status" %in% names(result))
  expect_equal(nrow(result), 2)
})

test_that("dilution_scan() skips genes with FDR <= threshold", {
  bulk_deg <- data.frame(
    gene = c("GENE1", "GENE2"),
    beta = c(-0.084, 0.020),
    se = c(0.063, 0.050),
    fdr = c(0.450, 0.01)
  )
  sc_zscores <- c(GENE1 = 2.077, GENE2 = 1.5)

  result <- dilution_scan(
    bulk_deg = bulk_deg,
    sc_zscores = sc_zscores,
    cluster_fraction = 0.0695,
    condition_fraction = 0.631,
    n_cluster = 1332,
    null_fdr_threshold = 0.05
  )

  expect_equal(result$gate_status[1], "PASSED")
  expect_equal(result$gate_status[2], "FAILED")
})

test_that("dilution_scan() handles missing single-cell z-scores", {
  bulk_deg <- data.frame(
    gene = c("GENE1", "GENE2"),
    beta = c(-0.084, 0.020),
    se = c(0.063, 0.050),
    fdr = c(0.450, 0.200)
  )
  sc_zscores <- c(GENE1 = 2.077)  # GENE2 missing

  result <- dilution_scan(
    bulk_deg = bulk_deg,
    sc_zscores = sc_zscores,
    cluster_fraction = 0.0695,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_equal(result$gate_status[2], "SKIPPED")
  expect_true(is.na(result$sc_zscore[2]))
})

test_that("dilution_scan() handles zero z-scores", {
  bulk_deg <- data.frame(
    gene = c("GENE1", "GENE2"),
    beta = c(-0.084, 0.020),
    se = c(0.063, 0.050),
    fdr = c(0.450, 0.200)
  )
  sc_zscores <- c(GENE1 = 2.077, GENE2 = 0)  # Zero z-score

  result <- dilution_scan(
    bulk_deg = bulk_deg,
    sc_zscores = sc_zscores,
    cluster_fraction = 0.0695,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_equal(result$gate_status[2], "SKIPPED")
  expect_match(result$verdict[2], "zero|Zero")
})

test_that("dilution_scan() computes correct DDS and DI for passing genes", {
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

  expect_equal(result$gate_status[1], "PASSED")
  expect_true(result$dds[1] >= 0 && result$dds[1] <= 1)
  expect_true(result$di[1] >= 0 && result$di[1] <= 1)
  expect_true(nchar(result$verdict[1]) > 0)
})

test_that("dilution_scan() supports alternative column names", {
  bulk_deg <- data.frame(
    gene = c("GENE1"),
    log2fc = c(-0.084),
    se_log2fc = c(0.063),
    padj = c(0.450)
  )
  sc_zscores <- c(GENE1 = 2.077)

  result <- dilution_scan(
    bulk_deg = bulk_deg,
    sc_zscores = sc_zscores,
    cluster_fraction = 0.0695,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_equal(nrow(result), 1)
  expect_equal(result$gate_status, "PASSED")
})

test_that("dilution_scan() produces verdicts for all processed genes", {
  bulk_deg <- data.frame(
    gene = c("GENE1", "GENE2", "GENE3"),
    beta = c(-0.084, 0.020, 0.150),
    se = c(0.063, 0.050, 0.100),
    fdr = c(0.450, 0.200, 0.001)
  )
  sc_zscores <- c(GENE1 = 2.077, GENE2 = 1.5, GENE3 = 0.5)

  result <- dilution_scan(
    bulk_deg = bulk_deg,
    sc_zscores = sc_zscores,
    cluster_fraction = 0.0695,
    condition_fraction = 0.631,
    n_cluster = 1332
  )

  expect_true(!all(is.na(result$verdict)))
  expect_true(result$gate_status[3] == "FAILED")
  expect_true(result$gate_status[1] %in% c("PASSED", "UNDERPOWERED", "AMBIGUOUS"))
})
