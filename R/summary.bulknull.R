#' Summary method for bulknull objects
#'
#' Displays a concise one-line summary of the diagnosis.
#'
#' @param object A bulknull object
#' @param ... Additional arguments (ignored)
#' @export
summary.bulknull <- function(object, ...) {
  cat("\n=== BULKNULL SUMMARY ===\n\n")
  cat("Status:    ", object$gate_status, "\n")
  if (!is.null(object$dds) && is.numeric(object$dds)) {
    cat("DDS:       ", round(object$dds, 4), " (", object$dds_interpretation, ")\n")
  }
  cat("DI:        ", round(object$di, 4), " (", object$di_interpretation, ")\n")
  cat("Verdict:   ", object$verdict, "\n\n")
  invisible(object)
}
