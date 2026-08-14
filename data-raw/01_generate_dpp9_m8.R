# Generate dpp9_m8.rda from config/dilution_parameters.csv
# Ensures data object is always synchronized with parameter specifications
# DO NOT edit dpp9_m8.rda by hand

# Read parameters
config_file <- "config/dilution_parameters.csv"
if (!file.exists(config_file)) {
  stop("config/dilution_parameters.csv not found. Run from project root.")
}

params <- read.csv(config_file, stringsAsFactors = FALSE, colClasses = "character")

# Extract k=5 tissue-only case study values
dpp9_m8 <- list(
  bulk_beta = as.numeric(params[params$value == "-0.084210", "value"]),
  bulk_se = as.numeric(params[params$value == "0.062774", "value"]),
  bulk_fdr = as.numeric(params[params$value == "0.449745", "value"]),
  sc_zscore = as.numeric(params[params$value == "2.077409", "value"]),
  cluster_fraction = 1332 / 19175,  # n_cluster / n_total (exact ratio)
  condition_fraction = as.numeric(params[params$value == "0.631", "value"]),
  n_cluster = 1332,

  # Metadata
  gene = "DPP9",
  cluster = "M8",
  cluster_name = "Myeloid (DPP9/NLRP1/CASP1-high inflammatory macrophages)",
  cohort = "IPF lung tissue (k=5 tissue-only: 5 series, BAL excluded, GSE10667 excluded due to scaling artifacts)",
  description = "DPP9 expression in M8 myeloid subcluster from IPF lung tissue (k=5 tissue-only meta-analysis). Primary case study for dilution_score and diagnosability_index validation. Bulk null despite strong single-cell signal (z=2.077409) in small cluster (3.3% of total cells)."
)

# Verify structure
stopifnot(
  is.numeric(dpp9_m8$bulk_beta),
  is.numeric(dpp9_m8$bulk_se),
  abs(dpp9_m8$bulk_beta - (-0.084210)) < 1e-6,
  abs(dpp9_m8$bulk_se - 0.062774) < 1e-6,
  abs(dpp9_m8$sc_zscore - 2.077409) < 1e-6
)

# Save to package
usethis::use_data(dpp9_m8, overwrite = TRUE, internal = FALSE)

cat("Generated dpp9_m8.rda from config/dilution_parameters.csv\n")
cat("k=5 tissue-only (5 lung IPF series, GSE10667 excluded)\n")
cat("Source: config/dilution_parameters.csv\n")
cat("Verify: library(bulknull); data(dpp9_m8); dpp9_m8\n")
