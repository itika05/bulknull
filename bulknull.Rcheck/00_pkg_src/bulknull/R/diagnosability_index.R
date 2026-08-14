#' Diagnosability Index (DI)
#'
#' Compute statistical power to detect dilution under a one-sided or two-sided test.
#'
#' @param mu_dilution Numeric. Expected signal under dilution:
#'   mu = (cluster_frac * sc_zscore) / (bulk_se * sqrt(n_eff)).
#' @param alpha Numeric. Significance level (default 0.05).
#' @param sided Character. "one" for one-sided test (default, recommended for
#'   directional alternatives) or "two" for two-sided.
#'
#' @return List containing:
#'   - di: Diagnosability Index (power), in range 0 to 1
#'   - critical_value: Critical z-value used
#'   - alpha: Significance level
#'   - sided: Test type ("one" or "two")
#'   - formula: Human-readable formula string
#'
#' @details
#' DI estimates the probability of correctly detecting dilution if it truly
#' occurred at the specified effect size (mu_dilution).
#'
#' One-sided (default, recommended):
#'   - Critical value: z = qnorm(1 - alpha) = 1.644854 for alpha=0.05
#'   - Rationale: Dilution alternative is directional by construction
#'     (sign(mu_dilution) = sign(sc_zscore) since cluster_frac > 0)
#'   - Formula: DI = pnorm(mu - 1.644854)
#'
#' Two-sided (for reference):
#'   - Critical value: z = qnorm(1 - alpha/2) = 1.959964 for alpha=0.05
#'   - Formula: DI = pnorm(mu - 1.959964)
#'
#' Power categories:
#'   - DI < 0.3: Low power (undetectable at this effect size)
#'   - DI 0.3-0.8: Moderate power (detection possible but unreliable)
#'   - DI > 0.8: High power (reliable detection)
#'
#' @examples
#' # DPP9/M8 k=5 tissue-only
#' diagnosability_index(mu_dilution = 0.130536, alpha = 0.05, sided = "one")
#' # Returns: di = 0.0650 (very low power; dilution undetectable)
#'
#' @importFrom stats pnorm qnorm
#' @export
diagnosability_index <- function(mu_dilution, alpha = 0.05, sided = c("one", "two")) {
  sided <- match.arg(sided)

  # Compute critical value based on test type
  if (sided == "one") {
    z_critical <- qnorm(1 - alpha)  # One-sided: 1.644854 for alpha=0.05
  } else {
    z_critical <- qnorm(1 - alpha / 2)  # Two-sided: 1.959964 for alpha=0.05
  }

  # Compute DI as power under H1: z ~ N(mu, 1)
  di <- pnorm(mu_dilution - z_critical)

  # Return result with details for reproducibility
  return(list(
    di = di,
    critical_value = z_critical,
    alpha = alpha,
    sided = sided,
    formula = sprintf("Normal(mu - %.6f)", z_critical)
  ))
}
