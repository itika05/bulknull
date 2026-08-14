# Data object synchronization test
# Ensures dpp9_m8 data object reproduces frozen values from NUMBERS_TABLE.tsv

test_that("dpp9_m8 data object contains correct k=5 tissue-only values", {
  # Load the data object
  data(dpp9_m8, envir = environment())

  # Verify structure
  expect_true(is.list(dpp9_m8))
  expect_named(dpp9_m8, c("bulk_beta", "bulk_se", "bulk_fdr", "sc_zscore",
                           "cluster_fraction", "condition_fraction", "n_cluster",
                           "gene", "cluster", "cluster_name", "cohort", "description"))

  # Verify k=5 tissue-only values (to 1e-6 precision)
  expect_equal(dpp9_m8$bulk_beta, -0.084210, tolerance = 1e-6)
  expect_equal(dpp9_m8$bulk_se, 0.062774, tolerance = 1e-6)
  expect_equal(dpp9_m8$bulk_fdr, 0.449745, tolerance = 1e-6)
  expect_equal(dpp9_m8$sc_zscore, 2.077409, tolerance = 1e-6)
  expect_equal(dpp9_m8$cluster_fraction, 1332 / 19175, tolerance = 1e-6)
  expect_equal(dpp9_m8$condition_fraction, 0.631, tolerance = 1e-6)
  expect_equal(dpp9_m8$n_cluster, 1332)

  # Verify metadata
  expect_equal(dpp9_m8$gene, "DPP9")
  expect_equal(dpp9_m8$cluster, "M8")
  expect_true(grepl("IPF", dpp9_m8$cohort))
  expect_true(grepl("tissue-only", dpp9_m8$cohort))
})

test_that("dpp9_m8 reproduces frozen dilution_score values", {
  data(dpp9_m8, envir = environment())

  # Compute dilution_score using data object
  result <- dilution_score(
    bulk_beta = dpp9_m8$bulk_beta,
    bulk_se = dpp9_m8$bulk_se,
    bulk_fdr = dpp9_m8$bulk_fdr,
    sc_zscore = dpp9_m8$sc_zscore,
    cluster_fraction = dpp9_m8$cluster_fraction,
    condition_fraction = dpp9_m8$condition_fraction,
    n_cluster = dpp9_m8$n_cluster
  )

  # Verify frozen values (z-score 2.077409, tolerance 1e-5 for floating point)
  expect_equal(result$mu_dilution, 0.1305362, tolerance = 1e-5)
  expect_equal(result$z_bulk, -1.3414790, tolerance = 1e-5)
  expect_equal(result$dds_score, 0.4542207, tolerance = 1e-5)
})

test_that("dpp9_m8 reproduces frozen diagnosability_index values", {
  data(dpp9_m8, envir = environment())

  # Compute diagnosability_index using frozen mu_dilution
  di_one <- diagnosability_index(mu_dilution = 0.130536, alpha = 0.05, sided = "one")
  di_two <- diagnosability_index(mu_dilution = 0.130536, alpha = 0.05, sided = "two")

  # Verify frozen DI values (to 1e-5 precision)
  expect_equal(di_one$di, 0.064973, tolerance = 1e-5)
  expect_equal(di_two$di, 0.033668, tolerance = 1e-5)
})
