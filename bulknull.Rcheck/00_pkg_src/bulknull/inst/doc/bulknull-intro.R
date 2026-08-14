## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(bulknull)
data(dpp9_m8)

## -----------------------------------------------------------------------------
str(dpp9_m8)

## -----------------------------------------------------------------------------
result <- dilution_score(
  bulk_beta = dpp9_m8$bulk_beta,
  bulk_se = dpp9_m8$bulk_se,
  bulk_fdr = dpp9_m8$bulk_fdr,
  sc_zscore = dpp9_m8$sc_zscore,
  cluster_fraction = dpp9_m8$cluster_fraction,
  condition_fraction = dpp9_m8$condition_fraction,
  n_cluster = dpp9_m8$n_cluster
)

## -----------------------------------------------------------------------------
cat("DDS Score:", result$dds_score, "\n")
cat("Interpretation:", result$summary, "\n")
cat("Applicable:", result$applicable, "\n")
cat("Note:", result$applicability_note, "\n")

## -----------------------------------------------------------------------------
# This case is applicable (bulk FDR > 0.05)
if (!result$applicable) {
  warning("This case is not applicable for DDS analysis.")
  warning(result$applicability_note)
}

## -----------------------------------------------------------------------------
# Compute mu_dilution from the components (same formula as inside dilution_score)
n_eff_sc <- dpp9_m8$condition_fraction * (1 - dpp9_m8$condition_fraction) * dpp9_m8$n_cluster
d_sc <- dpp9_m8$sc_zscore / sqrt(n_eff_sc)
d_bulk_expected <- dpp9_m8$cluster_fraction * d_sc
mu_dilution <- d_bulk_expected / dpp9_m8$bulk_se

# Compute DI
di_result <- diagnosability_index(
  mu_dilution = mu_dilution,
  alpha = 0.05,
  sided = "one"
)

## -----------------------------------------------------------------------------
cat("Diagnosability Index:", di_result$di, "\n")
if (di_result$di < 0.3) {
  cat("Low power: Study design may be inadequate for dilution detection.\n")
} else if (di_result$di > 0.8) {
  cat("High power: Study design provides good power for dilution detection.\n")
} else {
  cat("Moderate power: Study design has moderate power for dilution detection.\n")
}

## -----------------------------------------------------------------------------
# Example: Create a simple table of inputs
input_data <- data.frame(
  gene = c("DPP9", "SPP1", "CXCL10"),
  bulk_beta = c(-0.030, 0.150, 0.200),
  bulk_se = rep(0.060, 3),
  bulk_fdr = c(0.821, 0.005, 0.001),
  sc_zscore = c(2.077, 3.245, 2.890),
  cluster_fraction = c(0.0695, 0.1200, 0.0850),
  condition_fraction = rep(0.631, 3),
  n_cluster = c(1332, 2100, 1800)
)

# Apply DDS to each row (in practice, use sapply or mapply)
# results <- apply(input_data, 1, function(row) {
#   dilution_score(
#     bulk_beta = as.numeric(row["bulk_beta"]),
#     bulk_se = as.numeric(row["bulk_se"]),
#     bulk_fdr = as.numeric(row["bulk_fdr"]),
#     sc_zscore = as.numeric(row["sc_zscore"]),
#     cluster_fraction = as.numeric(row["cluster_fraction"]),
#     condition_fraction = as.numeric(row["condition_fraction"]),
#     n_cluster = as.numeric(row["n_cluster"])
#   )
# })

## -----------------------------------------------------------------------------
sessionInfo()

