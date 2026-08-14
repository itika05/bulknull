#' One-Call Wrapper for Bulk-Null Dilution Diagnostics
#'
#' Unified interface that runs the applicability gate, computes the Dilution
#' Discordance Score (DDS), and the Diagnosability Index (DI) in a single call.
#' Returns a formatted S3 object summarizing the complete analysis.
#'
#' @param bulk_beta numeric; effect size estimate from bulk analysis
#' @param bulk_se numeric; standard error of bulk_beta
#' @param bulk_fdr numeric; FDR (or p-value) from bulk test
#' @param d_sc numeric; standardized effect size in single-cell cluster (primary input)
#' @param sc_zscore numeric; z-score of signal in target cluster (backward compatibility)
#' @param cluster_fraction numeric; fraction of total cells in target cluster
#' @param condition_fraction numeric; fraction of cells in case condition
#' @param r numeric; relative expression (mean in condition / mean across compartment).
#'            Default 1 (no relative expression adjustment).
#' @param n_cluster numeric; number of cells in target cluster (required if sc_zscore supplied)
#' @param n_eff_basis character; basis for sc_zscore conversion: "cell" or "donor"
#' @param n_donors numeric; number of donors (required if n_eff_basis = "donor")
#' @param null_fdr_threshold numeric; FDR threshold for applicability gate (default 0.05)
#' @param on_violation character; "error", "warn", or "pass" (default "error")
#' @param alpha numeric; significance level for DI (default 0.05)
#' @param sided character; "one" or "two" for DI computation (default "one")
#' @param verbose logical; if TRUE, include type_s, type_m, and detailed DDS. Default FALSE.
#'
#' @return An S3 object of class "bulknull" containing:
#'   \item{mu_dilution}{Expected bulk z-score under dilution hypothesis}
#'   \item{z_bulk}{Observed bulk z-score}
#'   \item{w}{Mixture weight used}
#'   \item{di}{Diagnosability Index value}
#'   \item{r}{Relative expression parameter used}
#'   \item{verdict}{Interpretive verdict}
#'   \item{inputs}{List of all input parameters}
#'   \item{applicability}{List with $applicable and $note}
#'   \item{type_s}{Type S error rate (when verbose=TRUE)}
#'   \item{type_m}{Type M error ratio (when verbose=TRUE)}
#'   \item{dds}{DDS value (when verbose=TRUE only)}
#'
#' @details
#' The bulknull workflow is:
#' 1. Check applicability: bulk FDR > threshold (default 0.05)
#' 2. If applicable, compute DDS (posterior probability of dilution)
#' 3. Compute DI (power to detect dilution)
#' 4. Optionally compute type_s and type_m under model uncertainty
#'
#' The mixture weight w = (f*r)/(f*r + (1-f)) incorporates relative expression.
#' At r=1, this reduces to w=f, recovering the standard dilution model.
#'
#' @examples
#' # Example: DPP9 in M8 (direct d_sc input, recommended)
#' result <- bulknull(
#'   bulk_beta = -0.0297,
#'   bulk_se = 0.0595,
#'   bulk_fdr = 0.821,
#'   d_sc = 0.117962,
#'   cluster_fraction = 0.0695,
#'   condition_fraction = 0.631,
#'   alpha = 0.05,
#'   sided = "one"
#' )
#'
#' print(result)
#'
#' @importFrom stats qnorm pnorm dnorm
#' @export
bulknull <- function(bulk_beta, bulk_se, bulk_fdr, d_sc = NULL, sc_zscore = NULL,
                      cluster_fraction, condition_fraction, r = 1,
                      n_cluster = NULL,
                      n_eff_basis = c("cell", "donor"), n_donors = NULL,
                      null_fdr_threshold = 0.05,
                      on_violation = c("error", "warn", "pass"),
                      alpha = 0.05, sided = c("one", "two"), verbose = FALSE) {

  on_violation <- match.arg(on_violation)
  sided <- match.arg(sided)
  n_eff_basis <- match.arg(n_eff_basis)

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
    d_sc = d_sc,
    sc_zscore = sc_zscore,
    cluster_fraction = cluster_fraction,
    condition_fraction = condition_fraction,
    r = r,
    n_cluster = n_cluster,
    n_eff_basis = n_eff_basis,
    n_donors = n_donors,
    null_fdr_threshold = null_fdr_threshold,
    on_violation = on_violation
  )

  # Step 3: Compute DI
  di_result <- diagnosability_index(
    mu_dilution = dds_result$mu_dilution,
    alpha = alpha,
    sided = sided
  )

  # Step 4: Compute type S and type M (under N(mu, 1))
  mu <- dds_result$mu_dilution
  c <- if (sided == "one") qnorm(1 - alpha) else qnorm(1 - alpha / 2)

  # P(significant) using closed form
  p_sig <- pnorm(-c - mu) + (1 - pnorm(c - mu))

  # Type S = P(sign error | significant)
  type_s <- if (p_sig > 0) pnorm(-c - mu) / p_sig else NA_real_

  # Type M = E[|z| | |z| > c] / mu using exact closed form
  # E[|z| | |z| > c] * P(sig) = mu*(1-Phi(c-mu)) + phi(c-mu) - mu*Phi(-c-mu) + phi(-c-mu)
  if (p_sig > 0 && mu > 0) {
    numerator <- mu * (1 - pnorm(c - mu)) + dnorm(c - mu) -
                 mu * pnorm(-c - mu) + dnorm(-c - mu)
    e_abs_sig <- numerator / p_sig
    type_m <- e_abs_sig / mu
  } else {
    e_abs_sig <- NA_real_
    type_m <- NA_real_
  }

  # Build verdict (simplified, no interpretation bands per TASK 4)
  verdict <- if (!applicability$applicable) {
    "FRAMEWORK NOT APPLICABLE: Bulk effect is significant (FDR <= threshold)"
  } else if (di_result$di > 0.8) {
    sprintf("ADEQUATE POWER (DI=%.3f): Can detect dilution if present", di_result$di)
  } else if (di_result$di > 0.3) {
    sprintf("MODERATE POWER (DI=%.3f): Uncertain power to detect dilution", di_result$di)
  } else {
    sprintf("LOW POWER (DI=%.3f): Underpowered to detect dilution", di_result$di)
  }

  # Compute interpretations (needed by downstream functions)
  dds_interpretation <- dds_result$summary

  di_interpretation <- if (di_result$di > 0.8) {
    "High power"
  } else if (di_result$di > 0.3) {
    "Moderate power"
  } else {
    "Low power"
  }

  # Build output list (DDS removed from default per TASK 3)
  result <- list(
    inputs = list(
      bulk_beta = bulk_beta,
      bulk_se = bulk_se,
      bulk_fdr = bulk_fdr,
      d_sc = dds_result$d_sc,
      sc_zscore = sc_zscore,
      cluster_fraction = cluster_fraction,
      condition_fraction = condition_fraction,
      r = r,
      n_cluster = n_cluster,
      n_eff_basis = n_eff_basis,
      n_donors = n_donors,
      alpha = alpha,
      sided = sided
    ),
    applicability = list(
      applicable = applicability$applicable,
      note = applicability$applicability_note
    ),
    mu_dilution = dds_result$mu_dilution,
    z_bulk = dds_result$z_bulk,
    w = dds_result$w,
    d_sc = dds_result$d_sc,
    dds_interpretation = dds_interpretation,
    di = di_result$di,
    di_interpretation = di_interpretation,
    r = r,
    verdict = verdict,
    gate_status = if (applicability$applicable) "PASSED" else "FAILED",
    # Store full objects for inspection
    .dds_full = dds_result,
    .di_full = di_result
  )

  # Add type S, type M, DDS, and intermediate values only when verbose=TRUE
  if (verbose) {
    result$type_s <- type_s
    result$type_m <- type_m
    result$e_abs_sig <- e_abs_sig
    result$p_sig <- p_sig
    result$dds <- dds_result$dds_score
  }

  structure(result, class = "bulknull")
}
