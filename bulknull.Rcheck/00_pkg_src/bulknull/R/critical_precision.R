#' Critical Precision for Target Diagnosability Index
#'
#' Compute the minimum bulk measurement precision (SE) required to achieve a target
#' diagnosability index at any cluster fraction <= 1. A target DI is attainable if and
#' only if bulk_se <= s_critical, which depends only on d_sc, target_di, and alpha.
#'
#' The relationship is independent of the pooled cohort size N; it depends only on
#' the effect size standardization scale and the statistical threshold.
#'
#' @param d_sc numeric; standardized effect size in single-cell cluster
#' @param bulk_se numeric; observed bulk measurement standard error
#' @param target_di numeric; target diagnosability index (0 to 1, default 0.8)
#' @param alpha numeric; significance level (default 0.05)
#' @param sided character; "one" or "two" for one-sided or two-sided test (default "one")
#'
#' @return A list with elements:
#'   \item{mu_required}{Required effect size mu to achieve target DI}
#'   \item{s_critical}{Critical bulk SE; target DI attainable iff bulk_se <= s_critical}
#'   \item{attainable}{Logical; TRUE if target DI is attainable at current bulk_se}
#'   \item{min_detectable_d}{Smallest within-cluster effect reachable at f=1.0}
#'   \item{f_required_at_d}{Function f_required_at_d(d) returning required cluster fraction}
#'   \item{call}{The matched function call}
#'
#' @details
#' The condition for attainability is:
#'   bulk_se <= d_sc / (z_critical + qnorm(target_di))
#'
#' This formula is independent of N (the cohort size): it depends only on the
#' statistical threshold and effect-size scale.
#'
#' The minimum detectable effect at f=1 is:
#'   min_d = mu_required * bulk_se
#'
#' For a hypothesised within-cluster effect d, the required cluster fraction is:
#'   f_required = mu_required * bulk_se / d
#'
#' @examples
#' result <- critical_precision(
#'   d_sc = 0.117962,
#'   bulk_se = 0.062774,
#'   target_di = 0.8,
#'   alpha = 0.05,
#'   sided = "one"
#' )
#'
#' print(result$attainable)  # FALSE: bulk_se exceeds s_critical
#' print(result$min_detectable_d)
#'
#' @importFrom stats qnorm
#' @export
critical_precision <- function(d_sc, bulk_se, target_di = 0.8, alpha = 0.05,
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

  # Compute critical SE
  s_critical <- d_sc / mu_required

  # Attainability check
  attainable <- bulk_se <= s_critical

  # Min detectable effect at f=1
  min_detectable_d <- mu_required * bulk_se

  # Function to compute required cluster fraction for a given d
  f_required_at_d <- function(d = 0.8) {
    if (d == 0) return(Inf)
    mu_required * bulk_se / d
  }

  structure(
    list(
      mu_required = mu_required,
      s_critical = s_critical,
      attainable = attainable,
      min_detectable_d = min_detectable_d,
      f_required_at_d = f_required_at_d,
      call = match.call()
    ),
    class = "critical_precision"
  )
}

#' Print method for critical_precision objects
#'
#' @param x A critical_precision object
#' @param ... Additional arguments (ignored)
#' @export
print.critical_precision <- function(x, ...) {
  cat("=== Critical Precision Analysis ===\n\n")
  cat("Target Diagnosability Index: ", x$call$target_di, "\n")
  cat("Required mu: ", x$mu_required, "\n")
  cat("Critical bulk SE: ", x$s_critical, "\n")
  cat("Attainable: ", x$attainable, "\n")
  cat("Min detectable effect (at f=1): ", x$min_detectable_d, "\n")
  invisible(x)
}
