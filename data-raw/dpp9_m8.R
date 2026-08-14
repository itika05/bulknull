# Data preparation script for dpp9_m8 dataset
# This script creates the example dataset used in vignettes and tests

dpp9_m8 <- list(
  bulk_beta = -0.0297,
  bulk_se = 0.0595,
  bulk_fdr = 0.821,
  sc_zscore = 2.077,
  cluster_fraction = 1332 / 19175,
  condition_fraction = 0.631,
  n_cluster = 1332,
  gene = "DPP9",
  cluster = "M8",
  cluster_name = "Inflammatory Macrophages",
  cohort = "IPF",
  description = "DPP9 expression in M8 (inflammatory macrophages) subcluster from IPF lung tissue"
)

usethis::use_data(dpp9_m8, overwrite = TRUE)
