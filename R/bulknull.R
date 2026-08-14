#' One-Call Wrapper for Bulk-Null Dilution Diagnostics
#'
#' Unified interface that runs the applicability gate, computes the Dilution
#' Discordance Score (DDS), and the Diagnosability Index (DI) in a single call.
#' Returns a formatted S3 object summarizing the complete analysis.
#'
#' @param bulk_beta numeric; effect size estimate from bulk analysis
#' @param bulk_se numeric; standard error of bulk_beta
#' @param bulk_fdr numeric; FDR (or p-value) from bulk test
#' @param sc_zscore numeric; z-score of signal in target single-cell cluster
#' @param cluster_fraction numeric; fraction of total cells in target cluster
#' @param condition_fraction numeric; fraction of cells in case condition
#' @param n_cluster numeric; number of cells in target cluster
#' @param null_fdr_threshold numeric; FDR threshold for applicability gate (default 0.05)
#' @param on_violation character; "error", "warn", or "pass" (default "error")
#' @param alpha numeric; significance level for DI (default 0.05)
#' @param sided character; "one" or "two" for DI computation (default "one")
#'
#' @return An S3 object of class "bulknull" containing:
#'   \item{inputs}{List of input parameters}
#'   \item{applicability}{List with $applicable (logical) and $note (character)}
#'   \item{dds}{DDS (Dilution Discordance Score) value}
#'   \item{mu_dilution}{Expected bulk z-score under dilution hypothesis}
#'   \item{z_bulk}{Observed bulk z-score}
#'   \item{di}{Diagnosability Index value}
#'   \item{verdict}{Interpretive string (e.g., "Ambiguous")}
#'   \item{gate_status}{Character indicating whether input passed the gate}
#'
#' @details
#' The bulknull workflow is:
#' 1. Check applicability: bulk FDR > threshold (default 0.05)
#' 2. If applicable, compute DDS (posterior probability of dilution)
#' 3. Compute DI (power to detect dilution)
#' 4. Return formatted S3 object with interpretation bands
#'
#' The function preserves the full output of dilution_score() and
#' diagnosability_index() internally for inspection if needed.
#'
#' @examples
#' # Example: DPP9 in M8 myeloid cells from IPF cohort
#' result <- bulknull(
#'   bulk_beta = -0.084210,
#'   bulk_se = 0.062774,
#'   bulk_fdr = 0.449745,
#'   sc_zscore = 2.077409,
#'   cluster_fraction = 1332 / 19175,
#'   condition_fraction = 0.631,
#'   n_cluster = 1332,
#'   alpha = 0.05,
#'   sided = "one"
#' )
#'
#' print(result)
#'
#' @export
bulknull <- function(bulk_beta, bulk_se, bulk_fdr, sc_zscore, cluster_fraction,
                      condition_fraction, n_cluster, null_fdr_threshold = 0.05,
                      on_violation = c("error", "warn", "pass"),
                      alpha = 0.05, sided = c("one", "two")) {

  on_violation <- match.arg(on_violation)
  sided <- match.arg(sided)

  # Step 1: Check applicability
  applicability <- check_dds_applicability(
    bulk_fdr = bulk_fdr,
    null_fdr_threshold = null_fdr_threshold,
    on_violation = on_violation
  )

  # Step 2: Compute DDS
  dds_result <- dilution_score(
    bulk_beta = bulk_beta,
    bulk_se = bulk_se,
    bulk_fdr = bulk_fdr,
    sc_zscore = sc_zscore,
    cluster_fraction = cluster_fraction,
    condition_fraction = condition_fraction,
    n_cluster = n_cluster,
    null_fdr_threshold = null_fdr_threshold,
    on_violation = on_violation
  )

  # Step 3: Compute DI
  di_result <- diagnosability_index(
    mu_dilution = dds_result$mu_dilution,
    alpha = alpha,
    sided = sided
  )

  # Step 4: Interpret DDS and DI
  dds_interpretation <- if (!applicability$applicable) {
    "NOT APPLICABLE"
  } else if (dds_result$dds_score > 0.7) {
    "Strong evidence for dilution"
  } else if (dds_result$dds_score > 0.5) {
    "Weak-to-moderate evidence for dilution"
  } else if (dds_result$dds_score > 0.3) {
    "Weak evidence for dilution"
  } else {
    "Minimal evidence for dilution; consistent with genuine bulk null"
  }

  di_interpretation <- if (di_result$di < 0.3) {
    "Low power (underpowered)"
  } else if (di_result$di < 0.8) {
    "Moderate power (uncertain)"
  } else {
    "High power (well-powered)"
  }

  # Build verdict
  verdict <- if (!applicability$applicable) {
    "FRAMEWORK NOT APPLICABLE: Bulk effect is significant (FDR <= threshold)"
  } else if (dds_result$dds_score > 0.7 && di_result$di > 0.3) {
    "LIKELY DILUTION: Strong DDS evidence and adequate power"
  } else if (dds_result$dds_score > 0.5 && di_result$di > 0.3) {
    "POSSIBLE DILUTION: Moderate DDS evidence with adequate power"
  } else if (di_result$di < 0.3) {
    "UNDERPOWERED: Study lacks power to detect dilution; DDS result inconclusive"
  } else {
    "AMBIGUOUS: Weak DDS evidence; consistent with genuine null or underpowered dilution"
  }

  # Return S3 object
  structure(
    list(
      inputs = list(
        bulk_beta = bulk_beta,
        bulk_se = bulk_se,
        bulk_fdr = bulk_fdr,
        sc_zscore = sc_zscore,
        cluster_fraction = cluster_fraction,
        condition_fraction = condition_fraction,
        n_cluster = n_cluster
      ),
      applicability = list(
        applicable = applicability$applicable,
        note = applicability$applicability_note
      ),
      dds = dds_result$dds_score,
      mu_dilution = dds_result$mu_dilution,
      z_bulk = dds_result$z_bulk,
      dds_interpretation = dds_interpretation,
      di = di_result$di,
      di_interpretation = di_interpretation,
      verdict = verdict,
      gate_status = if (applicability$applicable) "PASSED" else "FAILED",
      # Store full objects for inspection
      .dds_full = dds_result,
      .di_full = di_result
    ),
    class = "bulknull"
  )
}
