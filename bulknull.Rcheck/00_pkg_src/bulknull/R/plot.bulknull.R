#' Plot DI Landscape for Study Design
#'
#' Visualizes the Diagnosability Index as a function of cluster fraction and
#' bulk measurement precision (SE), marking the current study's position. Base
#' graphics only.
#'
#' @param x A bulknull object
#' @param sc_zscore numeric; single-cell z-score (same as input to bulknull)
#' @param condition_fraction numeric; case fraction in cohort
#' @param n_cluster numeric; cells in target cluster
#' @param alpha numeric; significance level (default 0.05)
#' @param sided character; "one" or "two" (default "one")
#' @param cluster_fraction_range numeric vector of length 2; range for y-axis
#' @param se_range numeric vector of length 2; range for x-axis (log scale)
#' @param grid_points numeric; number of grid points per dimension (default 50)
#' @param ... Additional arguments (ignored)
#'
#' @details
#' The landscape shows DI as a 2D function of cluster_fraction and bulk_se.
#' Higher DI (>0.8) indicates good power; lower DI (<0.3) indicates underpowered.
#' Contours mark DI = 0.3 and DI = 0.8 boundaries. The study is marked with a
#' red dot at (bulk_se, cluster_fraction) from the input bulknull object.
#'
#' @importFrom graphics image contour points legend
#' @importFrom grDevices heat.colors
#' @importFrom stats pnorm qnorm
#' @export
plot.bulknull <- function(x, sc_zscore = NULL, condition_fraction = NULL,
                          n_cluster = NULL, alpha = 0.05,
                          sided = c("one", "two"),
                          cluster_fraction_range = c(0.001, 1),
                          se_range = c(0.001, 0.5),
                          grid_points = 50,
                          ...) {

  sided <- match.arg(sided)

  # Extract from bulknull object if not provided
  if (is.null(sc_zscore)) sc_zscore <- x$inputs$sc_zscore
  if (is.null(condition_fraction)) condition_fraction <- x$inputs$condition_fraction
  if (is.null(n_cluster)) n_cluster <- x$inputs$n_cluster

  # Compute z_critical
  z_critical <- if (sided == "one") stats::qnorm(1 - alpha) else stats::qnorm(1 - alpha / 2)

  # Create grid
  cf_seq <- seq(cluster_fraction_range[1], cluster_fraction_range[2], length.out = grid_points)
  se_seq <- 10^seq(log10(se_range[1]), log10(se_range[2]), length.out = grid_points)

  # Compute DI on grid
  n_eff_sc <- condition_fraction * (1 - condition_fraction) * n_cluster
  d_sc <- sc_zscore / sqrt(n_eff_sc)

  di_matrix <- matrix(NA, nrow = length(se_seq), ncol = length(cf_seq))

  for (i in seq_along(se_seq)) {
    for (j in seq_along(cf_seq)) {
      d_bulk <- cf_seq[j] * d_sc
      mu_dil <- d_bulk / se_seq[i]
      di_matrix[i, j] <- stats::pnorm(mu_dil - z_critical)
    }
  }

  # Create plot
  graphics::image(se_seq, cf_seq, di_matrix,
                  xlab = "Bulk SE (log scale)", ylab = "Cluster fraction",
                  main = "Diagnosability Index Landscape",
                  log = "x", col = grDevices::heat.colors(100))

  # Add contours at DI = 0.3 and 0.8
  graphics::contour(se_seq, cf_seq, di_matrix, levels = c(0.3, 0.8),
                    add = TRUE, drawlabels = TRUE, col = "black", lwd = 1.5)

  # Mark current study position
  graphics::points(x$inputs$bulk_se, x$inputs$cluster_fraction,
                   col = "red", pch = 16, cex = 1.5)
  graphics::legend("bottomright", pch = 16, col = "red",
                   legend = "Current study")

  invisible(x)
}
