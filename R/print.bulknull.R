#' Print method for bulknull objects
#'
#' Displays a formatted verdict block showing inputs, computed statistics,
#' interpretations, and the final verdict.
#'
#' @param x A bulknull object
#' @param ... Additional arguments (ignored)
#' @export
print.bulknull <- function(x, ...) {
  cat("\n=== BULKNULL DIAGNOSTIC RESULT ===\n\n")

  # Input parameters
  cat("INPUT PARAMETERS:\n")
  cat("  Bulk effect (beta):           ", x$inputs$bulk_beta, "\n")
  cat("  Bulk SE:                      ", x$inputs$bulk_se, "\n")
  cat("  Bulk FDR:                     ", x$inputs$bulk_fdr, "\n")
  cat("  Single-cell z-score:          ", x$inputs$sc_zscore, "\n")
  cat("  Cluster fraction:             ", x$inputs$cluster_fraction, "\n")
  cat("  Condition fraction:           ", x$inputs$condition_fraction, "\n")
  cat("  Cluster size (n):             ", x$inputs$n_cluster, "\n")

  # Applicability gate
  cat("\nAPPLICABILITY GATE:\n")
  cat("  Status:                       ", x$gate_status, "\n")
  if (!x$applicability$applicable) {
    cat("  Note:                         ", x$applicability$note, "\n")
  }

  # Computed statistics
  cat("\nCOMPUTED STATISTICS:\n")
  cat("  mu (dilution):                ", x$mu_dilution, "\n")
  cat("  z (bulk):                     ", x$z_bulk, "\n")

  # DDS interpretation
  cat("\nDILUTION DISCORDANCE SCORE (DDS):\n")
  cat("  Value:                        ", x$dds, "\n")
  cat("  Interpretation:               ", x$dds_interpretation, "\n")
  cat("  (>0.7 strong | >0.5 moderate | >0.3 weak | <=0.3 minimal)\n")

  # DI interpretation
  cat("\nDIAGNOSABILITY INDEX (DI):\n")
  cat("  Value:                        ", x$di, "\n")
  cat("  Interpretation:               ", x$di_interpretation, "\n")
  cat("  (<0.3 low power | <0.8 moderate | >=0.8 high power)\n")

  # Verdict
  cat("\nVERDICT:\n")
  cat("  ", x$verdict, "\n")

  cat("\n")
  invisible(x)
}
