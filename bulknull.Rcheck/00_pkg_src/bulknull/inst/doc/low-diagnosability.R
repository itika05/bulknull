## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bulknull)
data(dpp9_m8)

## ----di_baseline--------------------------------------------------------------
n_eff_sc <- dpp9_m8$condition_fraction * (1 - dpp9_m8$condition_fraction) * dpp9_m8$n_cluster
d_sc <- dpp9_m8$sc_zscore / sqrt(n_eff_sc)
d_bulk_expected <- dpp9_m8$cluster_fraction * d_sc
mu_dilution <- d_bulk_expected / dpp9_m8$bulk_se

di_result <- diagnosability_index(
  mu_dilution = mu_dilution,
  alpha = 0.05,
  sided = "one"
)

cat("DI at baseline:\n")
cat("  cluster_fraction =", dpp9_m8$cluster_fraction, "\n")
cat("  bulk_se =", dpp9_m8$bulk_se, "\n")
cat("  DI =", round(di_result$di, 4), "\n")

## ----option1------------------------------------------------------------------
design_cf <- required_design(
  target_di = 0.8,
  sc_zscore = dpp9_m8$sc_zscore,
  condition_fraction = dpp9_m8$condition_fraction,
  n_cluster = dpp9_m8$n_cluster,
  cluster_fraction = dpp9_m8$cluster_fraction,
  bulk_se = dpp9_m8$bulk_se
)

cat(design_cf$scenario_fixed_cf, "\n")
cat("Required cluster fraction: ",
    round(design_cf$required_cluster_fraction_raw, 4),
    " (impossible; maximum = 1.0)\n", sep = "")

## ----option2------------------------------------------------------------------
design_di08 <- required_design(
  target_di = 0.8,
  sc_zscore = dpp9_m8$sc_zscore,
  condition_fraction = dpp9_m8$condition_fraction,
  n_cluster = dpp9_m8$n_cluster,
  cluster_fraction = dpp9_m8$cluster_fraction,
  bulk_se = dpp9_m8$bulk_se
)

cat(design_di08$scenario_fixed_se, "\n")
fold_improvement_08 <- dpp9_m8$bulk_se / design_di08$required_bulk_se
cat("Fold improvement in precision (SE reduction) needed: ",
    round(fold_improvement_08, 1), "x\n", sep = "")

## ----option3------------------------------------------------------------------
design_di05 <- required_design(
  target_di = 0.5,
  sc_zscore = dpp9_m8$sc_zscore,
  condition_fraction = dpp9_m8$condition_fraction,
  n_cluster = dpp9_m8$n_cluster,
  cluster_fraction = dpp9_m8$cluster_fraction,
  bulk_se = dpp9_m8$bulk_se
)

cat(design_di05$scenario_fixed_se, "\n")
fold_improvement_05 <- dpp9_m8$bulk_se / design_di05$required_bulk_se
cat("Fold improvement in precision (SE reduction) needed: ",
    round(fold_improvement_05, 1), "x\n", sep = "")

## -----------------------------------------------------------------------------
sessionInfo()

