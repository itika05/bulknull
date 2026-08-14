#' DDS Bounds Over Bulk Z-Score Range
#'
#' Compute the attainable interval of Dilution Discordance Scores (DDS) for a given
#' mu_dilution (effect size) as the bulk z-score z_bulk varies over a specified range.
#'
#' The DDS is a function of both mu_dilution and z_bulk. For fixed mu_dilution,
#' extreme values of z_bulk produce extreme DDS values. This function returns the
#' interval of DDS that could arise if the true bulk z-score were anywhere in the
#' specified range (typically between negative 1.96 and positive 1.96 for typical bulk analyses).
#'
#' Identity: DDS(z) = plogis(mu * z - mu^2/2)
#'   where plogis is the logistic CDF (inverse of logit).
#'
#' @param mu numeric; effect size mu under dilution hypothesis
#' @param bulk_se numeric; bulk measurement standard error (used only for min_detectable_d)
#' @param z_lo numeric; lower bound of bulk z-score range (default -1.96)
#' @param z_hi numeric; upper bound of bulk z-score range (default 1.96)
#'
#' @return A list with elements:
#'   \item{lower}{Lower bound of attainable DDS interval}
#'   \item{upper}{Upper bound of attainable DDS interval}
#'   \item{width}{Width of the interval (upper - lower)}
#'   \item{reaches_0.7}{Logical; TRUE if interval reaches DDS >= 0.7}
#'   \item{reaches_0.9}{Logical; TRUE if interval reaches DDS >= 0.9}
#'   \item{call}{The matched function call}
#'
#' @details
#' The DDS is computed as:
#'   DDS = plogis(mu * z - mu^2/2)
#'
#' For fixed mu, as z varies from z_lo to z_hi, the DDS varies accordingly.
#' The function returns the minimum and maximum attainable DDS values.
#'
#' A wide interval means many DDS values are possible (uncertainty in true z_bulk).
#' A narrow interval means few DDS values are possible (strong constraint on outcome).
#'
#' @examples
#' bounds <- dds_bounds(
#'   mu = 0.130536,
#'   bulk_se = 0.062774,
#'   z_lo = -1.96,
#'   z_hi = 1.96
#' )
#'
#' cat("DDS range:", bounds$lower, "to", bounds$upper, "\n")
#' cat("Interval width:", bounds$width, "\n")
#' cat("Reaches strong evidence (0.7)?", bounds$reaches_0.7, "\n")
#'
#' @importFrom stats plogis
#' @export
dds_bounds <- function(mu, bulk_se = NULL, z_lo = -1.96, z_hi = 1.96) {

  # Compute DDS at both bounds using the logistic CDF
  dds_lo <- stats::plogis(mu * z_lo - mu^2 / 2)
  dds_hi <- stats::plogis(mu * z_hi - mu^2 / 2)

  # Determine lower and upper of the interval
  lower <- min(dds_lo, dds_hi)
  upper <- max(dds_lo, dds_hi)

  # Width of interval
  width <- upper - lower

  # Check milestone evidence thresholds
  reaches_0.7 <- upper >= 0.7
  reaches_0.9 <- upper >= 0.9

  structure(
    list(
      lower = lower,
      upper = upper,
      width = width,
      reaches_0.7 = reaches_0.7,
      reaches_0.9 = reaches_0.9,
      call = match.call()
    ),
    class = "dds_bounds"
  )
}

#' Print method for dds_bounds objects
#'
#' @param x A dds_bounds object
#' @param ... Additional arguments (ignored)
#' @export
print.dds_bounds <- function(x, ...) {
  cat("=== DDS Bounds Analysis ===\n\n")
  cat("DDS interval: [", x$lower, ", ", x$upper, "]\n", sep = "")
  cat("Interval width: ", x$width, "\n", sep = "")
  cat("Reaches strong evidence (DDS >= 0.7)? ", x$reaches_0.7, "\n", sep = "")
  cat("Reaches very strong evidence (DDS >= 0.9)? ", x$reaches_0.9, "\n", sep = "")
  invisible(x)
}
