test_that("DDS increases with smaller cluster fraction (greater dilution)", {
  dds_large_cluster <- dilution_score(
    bulk_beta = 0.0,
    bulk_se = 0.06,
    bulk_fdr = 0.5,
    sc_zscore = 2.0,
    cluster_fraction = 0.20,
    condition_fraction = 0.5,
    n_cluster = 2000
  )$dds_score

  dds_small_cluster <- dilution_score(
    bulk_beta = 0.0,
    bulk_se = 0.06,
    bulk_fdr = 0.5,
    sc_zscore = 2.0,
    cluster_fraction = 0.02,
    condition_fraction = 0.5,
    n_cluster = 200
  )$dds_score

  expect_true(dds_small_cluster > dds_large_cluster)
})

test_that("DDS discriminates between weak and strong single-cell signals given bulk null", {
  dds_weak_sc <- dilution_score(
    bulk_beta = 0.0,
    bulk_se = 0.06,
    bulk_fdr = 0.5,
    sc_zscore = 1.0,
    cluster_fraction = 0.05,
    condition_fraction = 0.5,
    n_cluster = 500
  )$dds_score

  dds_strong_sc <- dilution_score(
    bulk_beta = 0.0,
    bulk_se = 0.06,
    bulk_fdr = 0.5,
    sc_zscore = 3.0,
    cluster_fraction = 0.05,
    condition_fraction = 0.5,
    n_cluster = 500
  )$dds_score

  expect_true(dds_strong_sc < dds_weak_sc)
})

test_that("DDS sensitivity changes with condition fraction (statistical power effect)", {
  dds_balanced <- dilution_score(
    bulk_beta = 0.0,
    bulk_se = 0.06,
    bulk_fdr = 0.5,
    sc_zscore = 2.0,
    cluster_fraction = 0.05,
    condition_fraction = 0.50,
    n_cluster = 500
  )$dds_score

  dds_imbalanced_high <- dilution_score(
    bulk_beta = 0.0,
    bulk_se = 0.06,
    bulk_fdr = 0.5,
    sc_zscore = 2.0,
    cluster_fraction = 0.05,
    condition_fraction = 0.90,
    n_cluster = 500
  )$dds_score

  dds_imbalanced_low <- dilution_score(
    bulk_beta = 0.0,
    bulk_se = 0.06,
    bulk_fdr = 0.5,
    sc_zscore = 2.0,
    cluster_fraction = 0.05,
    condition_fraction = 0.10,
    n_cluster = 500
  )$dds_score

  expect_true(dds_balanced > dds_imbalanced_high)
  expect_true(dds_balanced > dds_imbalanced_low)
})
