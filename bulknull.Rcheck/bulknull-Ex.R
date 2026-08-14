pkgname <- "bulknull"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "bulknull-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('bulknull')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("bulknull")
### * bulknull

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bulknull
### Title: One-Call Wrapper for Bulk-Null Dilution Diagnostics
### Aliases: bulknull

### ** Examples

# Example: DPP9 in M8 (direct d_sc input, recommended)
result <- bulknull(
  bulk_beta = -0.0297,
  bulk_se = 0.0595,
  bulk_fdr = 0.821,
  d_sc = 0.117962,
  cluster_fraction = 0.0695,
  condition_fraction = 0.631,
  alpha = 0.05,
  sided = "one"
)

print(result)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bulknull", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("check_dds_applicability")
### * check_dds_applicability

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: check_dds_applicability
### Title: Check DDS Applicability: Bulk Must Be Null
### Aliases: check_dds_applicability

### ** Examples

# Case 1: Bulk is non-significant (FDR = 0.15) - APPLICABLE
check_dds_applicability(bulk_fdr = 0.15, null_fdr_threshold = 0.05,
                        on_violation = "warn")

# Case 2: Bulk is significant (FDR = 0.001) - NOT APPLICABLE
check_dds_applicability(bulk_fdr = 0.001, null_fdr_threshold = 0.05,
                        on_violation = "warn")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("check_dds_applicability", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cohort_inflation")
### * cohort_inflation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cohort_inflation
### Title: Cohort Inflation Factor for Target Diagnosability Index
### Aliases: cohort_inflation

### ** Examples

inflation <- cohort_inflation(
  mu_observed = 0.130536,
  target_di = 0.8,
  alpha = 0.05,
  sided = "one"
)

cat("Required inflation:", inflation, "x\n")
cat("If current N=630, new N would be:", 630 * inflation, "\n")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cohort_inflation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("critical_precision")
### * critical_precision

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: critical_precision
### Title: Critical Precision for Target Diagnosability Index
### Aliases: critical_precision

### ** Examples

result <- critical_precision(
  d_sc = 0.117962,
  bulk_se = 0.062774,
  target_di = 0.8,
  alpha = 0.05,
  sided = "one"
)

print(result$attainable)  # FALSE: bulk_se exceeds s_critical
print(result$min_detectable_d)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("critical_precision", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dds_bounds")
### * dds_bounds

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dds_bounds
### Title: DDS Bounds Over Bulk Z-Score Range
### Aliases: dds_bounds

### ** Examples

bounds <- dds_bounds(
  mu = 0.130536,
  bulk_se = 0.062774,
  z_lo = -1.96,
  z_hi = 1.96
)

cat("DDS range:", bounds$lower, "to", bounds$upper, "\n")
cat("Interval width:", bounds$width, "\n")
cat("Reaches strong evidence (0.7)?", bounds$reaches_0.7, "\n")




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dds_bounds", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("diagnosability_index")
### * diagnosability_index

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: diagnosability_index
### Title: Diagnosability Index (DI)
### Aliases: diagnosability_index

### ** Examples

# DPP9/M8 k=5 tissue-only
diagnosability_index(mu_dilution = 0.130536, alpha = 0.05, sided = "one")
# Returns: di = 0.0650 (very low power; dilution undetectable)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("diagnosability_index", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dilution_scan")
### * dilution_scan

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dilution_scan
### Title: Vectorized Dilution Scan Across Multiple Genes
### Aliases: dilution_scan

### ** Examples

# Example: scan a small DEG table with mock z-scores
bulk_deg <- data.frame(
  gene = c("GENE1", "GENE2", "GENE3"),
  beta = c(-0.084, 0.020, 0.150),
  se = c(0.063, 0.050, 0.100),
  fdr = c(0.450, 0.05, 0.001)
)
sc_zscores <- c(GENE1 = 2.077, GENE2 = 1.5, GENE3 = 0.5)

result <- dilution_scan(
  bulk_deg = bulk_deg,
  sc_zscores = sc_zscores,
  cluster_fraction = 0.0695,
  condition_fraction = 0.631,
  n_cluster = 1332
)
print(result)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dilution_scan", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dilution_score")
### * dilution_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dilution_score
### Title: Dilution Discordance Score (DDS)
### Aliases: dilution_score

### ** Examples

# Direct d_sc input (recommended)
dilution_score(
  bulk_beta = -0.0297,
  bulk_se = 0.0595,
  bulk_fdr = 0.821,
  d_sc = 0.117962,
  cluster_fraction = 0.0695,
  condition_fraction = 0.631
)

# Backward compatible: sc_zscore with cell-level basis
dilution_score(
  bulk_beta = -0.0297,
  bulk_se = 0.0595,
  bulk_fdr = 0.821,
  sc_zscore = 2.077,
  cluster_fraction = 0.0695,
  condition_fraction = 0.631,
  n_cluster = 1332,
  n_eff_basis = "cell"
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dilution_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("dpp9_m8")
### * dpp9_m8

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: dpp9_m8
### Title: DPP9/M8 Example Dataset
### Aliases: dpp9_m8
### Keywords: datasets

### ** Examples

data(dpp9_m8)
str(dpp9_m8)

# Compute Dilution Discordance Score
result <- dilution_score(
  bulk_beta = dpp9_m8$bulk_beta,
  bulk_se = dpp9_m8$bulk_se,
  bulk_fdr = dpp9_m8$bulk_fdr,
  sc_zscore = dpp9_m8$sc_zscore,
  cluster_fraction = dpp9_m8$cluster_fraction,
  condition_fraction = dpp9_m8$condition_fraction,
  n_cluster = dpp9_m8$n_cluster
)
print(result$dds_score)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("dpp9_m8", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("required_design")
### * required_design

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: required_design
### Title: Required Study Design to Achieve Target Diagnosability Index
### Aliases: required_design

### ** Examples

# Example: k=5 tissue-only case study
# Current DI is only 0.065; what would we need for DI = 0.3?

design <- required_design(
  target_di = 0.3,
  sc_zscore = 2.077409,
  condition_fraction = 0.631,
  n_cluster = 1332,
  cluster_fraction = 1332 / 19175,
  bulk_se = 0.062774,
  alpha = 0.05,
  sided = "one"
)

print(design)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("required_design", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
