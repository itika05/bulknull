#' Required Study Design to Achieve Target Diagnosability Index
#'
#' Inverts the diagnosability index formula to find the study parameters needed
#' to achieve a target power (DI). Can solve for either cluster_fraction (at fixed
#' bulk_se) or bulk_se (at fixed cluster_fraction).
#'
#' @param target_di numeric; target diagnosability index (0 to 1)
#' @param d_sc numeric; standardized effect size in single-cell cluster (primary input)
#' @param sc_zscore numeric; single-cell z-score (backward compatibility)
#' @param condition_fraction numeric; fraction of cells in case condition (0 to 1)
#' @param n_cluster numeric; number of cells in target cluster (required if sc_zscore supplied)
#' @param n_eff_basis character; basis for conversion: "cell" or "donor"
#' @param n_donors numeric; number of donors (required if n_eff_basis = "donor")
#' @param cluster_fraction numeric; fraction of total cells in target cluster
#' @param bulk_se numeric; standard error of bulk effect estimate
#' @param alpha numeric; significance level (default 0.05)
#' @param sided character; "one" or "two" for one-sided or two-sided test (default "one")
#'
#' @return A list with class "required_design" containing:
#'   \item{target_di}{Requested diagnosability index}
#'   \item{achieved_di_at_fixed_se}{DI achieved if solving for cluster_fraction}
#'   \item{required_cluster_fraction}{Cluster fraction needed to reach target DI (at fixed bulk_se)}
#'   \item{required_bulk_se}{Bulk SE needed to reach target DI (at fixed cluster_fraction)}
#'   \item{achieved_di_at_fixed_cf}{DI achieved if solving for bulk_se}
#'   \item{scenario_fixed_se}{Message describing the fixed-SE scenario}
#'   \item{scenario_fixed_cf}{Message describing the fixed-CF scenario}
#'   \item{warnings}{Character vector of warnings if targets are unattainable}
#'
#' @details
#' The diagnosability index is defined as: DI = Normal_CDF(mu - z_critical)
#' where Normal_CDF is the standard normal CDF and z_critical depends on alpha and sided.
#'
#' This function inverts the formula in two ways:
#' 1. Solves for cluster_fraction given bulk_se, sc_zscore, and condition_fraction
#' 2. Solves for bulk_se given cluster_fraction, sc_zscore, and condition_fraction
#'
#' If the required cluster_fraction exceeds 1 (unattainable), the function returns NA
#' and includes a warning message.
#'
#' @examples
#' # Example: k=5 tissue-only case study
#' # Current DI is only 0.065; what would we need for DI = 0.3?
#'
#' design <- required_design(
#'   target_di = 0.3,
#'   sc_zscore = 2.077409,
#'   condition_fraction = 0.631,
#'   n_cluster = 1332,
#'   cluster_fraction = 1332 / 19175,
#'   bulk_se = 0.062774,
#'   alpha = 0.05,
#'   sided = "one"
#' )
#'
#' print(design)
#'
#' @export
required_design <- function(target_di, d_sc = NULL, sc_zscore = NULL,
                             condition_fraction, n_cluster = NULL,
                             n_eff_basis = c("cell", "donor"), n_donors = NULL,
                             cluster_fraction, bulk_se, alpha = 0.05,
                             sided = c("one", "two")) {

  sided <- match.arg(sided)
  n_eff_basis <- match.arg(n_eff_basis)

  # Handle d_sc vs sc_zscore: d_sc takes priority
  if (!is.null(d_sc) && !is.null(sc_zscore)) {
    warning("Both d_sc and sc_zscore supplied; using d_sc and ignoring sc_zscore")
  }

  if (is.null(d_sc)) {
    # Convert sc_zscore to d_sc
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

  # Compute z_critical
  if (sided == "one") {
    z_critical <- stats::qnorm(1 - alpha)
  } else {
    z_critical <- stats::qnorm(1 - alpha / 2)
  }

  # Solve for mu_dilution using inverse normal CDF
  # DI = Normal_CDF(mu - z_critical)
  # mu = Normal_CDF_inverse(DI) + z_critical = qnorm(DI) + z_critical
  mu_dilution_required <- stats::qnorm(target_di) + z_critical

  # Scenario 1: Solve for cluster_fraction at fixed bulk_se
  # μ = (cluster_frac * d_sc) / bulk_se
  # cluster_frac = μ * bulk_se / d_sc
  required_cf_raw <- mu_dilution_required * bulk_se / d_sc
  required_cf <- required_cf_raw

  # Check if cluster_fraction is attainable (≤ 1)
  warning_cf <- NULL
  if (required_cf > 1.0) {
    warning_cf <- sprintf("Target DI %.4f requires cluster_fraction %.4f > 1.0 (impossible).", target_di, required_cf_raw)
    required_cf <- NA
  }

  # Recompute achieved DI if CF was attainable
  if (!is.na(required_cf)) {
    d_bulk_expected <- required_cf * d_sc
    mu_dilution_achieved_cf <- d_bulk_expected / bulk_se
    achieved_di_cf <- stats::pnorm(mu_dilution_achieved_cf - z_critical)
  } else {
    achieved_di_cf <- NA
  }

  # Scenario 2: Solve for bulk_se at fixed cluster_fraction
  # We need: μ = (CF * d_sc) / bulk_se_required
  # bulk_se_required = (CF * d_sc) / μ
  d_bulk_expected <- cluster_fraction * d_sc
  required_se <- d_bulk_expected / mu_dilution_required

  # Recompute achieved DI with this SE
  achieved_di_se <- stats::pnorm(mu_dilution_required - z_critical)

  # Build result
  scenario_msg_se <- sprintf(
    "To achieve DI = %.4f at cluster_fraction = %.6f:\n  - Required bulk_se: %.6f\n  - Current bulk_se: %.6f",
    achieved_di_se, cluster_fraction, required_se, bulk_se
  )

  scenario_msg_cf <- if (!is.na(required_cf)) {
    sprintf(
      "To achieve DI = %.4f at bulk_se = %.6f:\n  - Required cluster_fraction: %.4f",
      achieved_di_cf, bulk_se, required_cf
    )
  } else {
    sprintf("Target DI = %.4f is unattainable (would require cluster_fraction > 1.0).", target_di)
  }

  warnings <- c()
  if (!is.null(warning_cf)) warnings <- c(warnings, warning_cf)

  structure(
    list(
      target_di = target_di,
      required_cluster_fraction = required_cf,
      required_cluster_fraction_raw = required_cf_raw,
      achieved_di_at_fixed_se = achieved_di_cf,
      required_bulk_se = required_se,
      achieved_di_at_fixed_cf = achieved_di_se,
      scenario_fixed_cf = scenario_msg_cf,
      scenario_fixed_se = scenario_msg_se,
      warnings = if (length(warnings) == 0) NULL else warnings
    ),
    class = "required_design"
  )
}

#' Print method for required_design objects
#'
#' @param x A required_design object
#' @param ... Additional arguments (ignored)
#' @export
print.required_design <- function(x, ...) {
  cat("=== Required Study Design for Target DI ===\n\n")
  cat("Target Diagnosability Index:", x$target_di, "\n\n")

  cat("SCENARIO 1: Solve for cluster_fraction (fixed bulk_se)\n")
  cat(x$scenario_fixed_cf, "\n\n")

  cat("SCENARIO 2: Solve for bulk_se (fixed cluster_fraction)\n")
  cat(x$scenario_fixed_se, "\n\n")

  if (!is.null(x$warnings)) {
    cat("WARNINGS:\n")
    for (w in x$warnings) {
      cat("  -", w, "\n")
    }
  }

  invisible(x)
}
