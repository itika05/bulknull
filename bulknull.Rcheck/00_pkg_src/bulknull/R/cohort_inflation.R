#' Cohort Inflation Factor for Target Diagnosability Index
#'
#' Compute the required inflation in cohort size (and thus reduction in bulk SE)
#' to achieve a target diagnosability index. The inflation factor is independent of
#' the pooled cohort size N and the bulk point estimate; it depends only on the
#' observed and required effect sizes mu.
#'
#' Identity: required inflation = (mu_required / mu_observed)^2
#'
#' @param mu_observed numeric; observed effect size (bulk z-score / sqrt(n_eff_bulk))
#' @param target_di numeric; target diagnosability index (0 to 1, default 0.8)
#' @param alpha numeric; significance level (default 0.05)
#' @param sided character; "one" or "two" for one-sided or two-sided test (default "one")
#'
#' @return A numeric scalar: the required inflation factor (>= 1.0 if target is attainable)
#'
#' @details
#' The relationship between bulk SE and cohort size is:
#'   SE ~ 1 / sqrt(N)
#'
#' To reduce SE by a factor k, cohort size must increase by k^2. Therefore,
#' if the current mu is mu_obs and we need mu_req, the inflation is:
#'   inflation = (mu_req / mu_obs)^2
#'
#' This computation is independent of:
#' - The absolute cohort size
#' - The bulk point estimate (only mu_obs matters)
#' - The cluster fraction or single-cell effect size
#'
#' An inflation of 50x means the cohort must grow 50 times in size.
#' An inflation of 50 with current N=630 implies new N ~ 31,500.
#'
#' @examples
#' inflation <- cohort_inflation(
#'   mu_observed = 0.130536,
#'   target_di = 0.8,
#'   alpha = 0.05,
#'   sided = "one"
#' )
#'
#' cat("Required inflation:", inflation, "x\n")
#' cat("If current N=630, new N would be:", 630 * inflation, "\n")
#'
#' @importFrom stats qnorm
#' @export
cohort_inflation <- function(mu_observed, target_di = 0.8, alpha = 0.05,
                              sided = c("one", "two")) {

  sided <- match.arg(sided)

  # Compute z_critical
  if (sided == "one") {
    z_critical <- stats::qnorm(1 - alpha)
  } else {
    z_critical <- stats::qnorm(1 - alpha / 2)
  }

  # Compute mu_required
  mu_required <- stats::qnorm(target_di) + z_critical

  # Inflation factor
  (mu_required / mu_observed)^2
}
