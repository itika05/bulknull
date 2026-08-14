#' Dilution Discordance Score (DDS)
#'
#' Compute the Bayesian posterior probability that an observed bulk-level effect
#' reflects true dilution of a single-cell signal.
#'
#' @param bulk_beta Numeric. Observed bulk effect size (Hedges' g or log2-fold-change).
#' @param bulk_se Numeric. Standard error of bulk effect size.
#' @param bulk_fdr Numeric. FDR-adjusted p-value from bulk analysis (or raw p-value).
#'                          Used for validation only; not directly in likelihood.
#' @param d_sc Numeric. Standardized effect size in the single-cell cluster (primary input).
#'                     When supplied, sc_zscore conversion is skipped.
#' @param sc_zscore Numeric. Single-cell z-score (backward compatibility). Ignored if d_sc supplied.
#'                           Must provide either d_sc or sc_zscore, not both unless testing.
#' @param cluster_fraction Numeric. Fraction of total tissue represented by the
#'                                  responsible subcluster (0 < f <= 1).
#' @param condition_fraction Numeric. Fraction of subcluster cells from the
#'                                    condition of interest (e.g., IPF, 0 <= p <= 1).
#' @param n_cluster Integer. Number of cells in the subcluster. Required if sc_zscore supplied.
#' @param n_eff_basis Character. Basis for sc_zscore conversion: "cell" or "donor".
#'                             Default "cell". Ignored if d_sc supplied.
#' @param n_donors Integer. Number of donors. Required if n_eff_basis = "donor".
#' @param null_fdr_threshold Numeric. FDR threshold for bulk to be considered "null".
#'                           Default: 0.05. Only used if bulk_fdr is provided.
#' @param on_violation Character. One of "error", "warn", "pass". Controls behavior
#'                     when bulk FDR is significant. Default: "warn".
#'
#' @return
#' A list with elements:
#'   - dds_score: Posterior probability of dilution (0 to 1).
#'   - z_bulk: Bulk z-score (beta / SE).
#'   - mu_dilution: Expected z under dilution model.
#'   - d_sc: Standardized effect size used (supplied or converted).
#'   - log_lik_h1: Log-likelihood under dilution model.
#'   - log_lik_h0: Log-likelihood under null model.
#'   - summary: Plain-text interpretation.
#'   - applicable: Logical. TRUE if case is applicable (bulk is null), FALSE if not.
#'   - applicability_note: Character string describing applicability status.
#'   - call: The matched function call.
#'
#' @examples
#' # Direct d_sc input (recommended)
#' dilution_score(
#'   bulk_beta = -0.0297,
#'   bulk_se = 0.0595,
#'   bulk_fdr = 0.821,
#'   d_sc = 0.117962,
#'   cluster_fraction = 0.0695,
#'   condition_fraction = 0.631
#' )
#'
#' # Backward compatible: sc_zscore with cell-level basis
#' dilution_score(
#'   bulk_beta = -0.0297,
#'   bulk_se = 0.0595,
#'   bulk_fdr = 0.821,
#'   sc_zscore = 2.077,
#'   cluster_fraction = 0.0695,
#'   condition_fraction = 0.631,
#'   n_cluster = 1332,
#'   n_eff_basis = "cell"
#' )
#'
#' @importFrom stats dnorm
#' @export
dilution_score <- function(
    bulk_beta,
    bulk_se,
    bulk_fdr,
    d_sc = NULL,
    sc_zscore = NULL,
    cluster_fraction,
    condition_fraction,
    n_cluster = NULL,
    n_eff_basis = c("cell", "donor"),
    n_donors = NULL,
    null_fdr_threshold = 0.05,
    on_violation = c("error", "warn", "pass")
) {

  on_violation <- match.arg(on_violation)
  n_eff_basis <- match.arg(n_eff_basis)

  # Input validation
  if (bulk_se <= 0) stop("bulk_se must be positive")
  if (!(0 < cluster_fraction && cluster_fraction <= 1)) {
    stop("cluster_fraction must be in (0, 1]")
  }
  if (!(0 <= condition_fraction && condition_fraction <= 1)) {
    stop("condition_fraction must be in [0, 1]")
  }

  # Handle d_sc vs sc_zscore: d_sc takes priority
  if (!is.null(d_sc) && !is.null(sc_zscore)) {
    warning("Both d_sc and sc_zscore supplied; using d_sc and ignoring sc_zscore")
  }

  if (is.null(d_sc)) {
    # Convert sc_zscore to d_sc using specified basis
    if (is.null(sc_zscore)) {
      stop("Must supply either d_sc or sc_zscore")
    }

    if (n_eff_basis == "cell") {
      if (is.null(n_cluster)) {
        stop("n_cluster required when n_eff_basis = 'cell'")
      }
      if (n_cluster < 1) stop("n_cluster must be >= 1")
      n_eff_sc <- condition_fraction * (1 - condition_fraction) * n_cluster
      d_sc <- sc_zscore / sqrt(n_eff_sc)
    } else if (n_eff_basis == "donor") {
      if (is.null(n_donors)) {
        stop("n_donors required when n_eff_basis = 'donor'")
      }
      if (n_donors < 2) stop("n_donors must be >= 2")
      n_eff_sc <- condition_fraction * (1 - condition_fraction) * n_donors
      d_sc <- sc_zscore / sqrt(n_eff_sc)
    }
  }

  # Applicability gate
  gate_result <- check_dds_applicability(
    bulk_fdr = bulk_fdr,
    null_fdr_threshold = null_fdr_threshold,
    on_violation = on_violation
  )

  # Compute bulk z-score and mu_dilution
  z_bulk <- bulk_beta / bulk_se
  d_bulk_expected <- cluster_fraction * d_sc
  mu_dilution <- d_bulk_expected / bulk_se

  # Likelihood ratio test
  log_lik_h1 <- dnorm(z_bulk, mean = mu_dilution, sd = 1, log = TRUE)
  log_lik_h0 <- dnorm(z_bulk, mean = 0, sd = 1, log = TRUE)

  # Convert to probability (posterior probability of H1)
  log_posterior_odds <- log_lik_h1 - log_lik_h0
  dds_score <- 1 / (1 + exp(-log_posterior_odds))

  # Interpretation
  if (dds_score > 0.9) {
    interpretation <- "STRONG evidence for dilution (dds > 0.9)"
  } else if (dds_score > 0.7) {
    interpretation <- "MODERATE-to-STRONG evidence for dilution (0.7 < dds <= 0.9)"
  } else if (dds_score > 0.5) {
    interpretation <- "WEAK-to-MODERATE evidence for dilution (0.5 < dds <= 0.7)"
  } else if (dds_score > 0.3) {
    interpretation <- "WEAK evidence for dilution (0.3 < dds <= 0.5)"
  } else {
    interpretation <- paste("MINIMAL evidence for dilution (dds <= 0.3);",
                            "consistent with genuine null")
  }

  # Return structured result
  return(
    list(
      dds_score = dds_score,
      z_bulk = z_bulk,
      mu_dilution = mu_dilution,
      d_sc = d_sc,
      log_lik_h1 = log_lik_h1,
      log_lik_h0 = log_lik_h0,
      summary = interpretation,
      applicable = gate_result$applicable,
      applicability_note = gate_result$applicability_note,
      call = match.call()
    )
  )
}
