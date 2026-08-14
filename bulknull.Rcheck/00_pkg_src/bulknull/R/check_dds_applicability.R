#' Check DDS Applicability: Bulk Must Be Null
#'
#' @description
#' Validates whether a case is applicable for DDS analysis.
#'
#' DDS answers the question: "Given a bulk-NULL result and strong single-cell signal,
#' is the null due to dilution or true absence of effect?"
#'
#' This function enforces the applicability boundary:
#' - **Applicable:** bulk FDR > null_fdr_threshold (non-significant bulk result)
#' - **Not Applicable:** bulk FDR <= null_fdr_threshold (significant bulk result)
#'
#' When not applicable, DDS is mathematically valid but answers the wrong question.
#'
#' @param bulk_fdr Numeric. FDR-adjusted p-value from bulk analysis.
#' @param null_fdr_threshold Numeric. FDR threshold for bulk to be considered "null".
#'                           Default: 0.05 (standard significance threshold).
#' @param on_violation Character. One of:
#'   - "error": Stop with informative error message (strict mode)
#'   - "warn": Issue warning but return results (permissive mode)
#'   - "pass": Return results silently (developer/testing mode)
#'
#' @return
#' A list with elements:
#'   - applicable: Logical. TRUE if bulk_fdr > null_fdr_threshold (NULL), FALSE otherwise.
#'   - applicability_note: Character string explaining the status.
#'   - violation_severity: Character. "none", "warning", or "error" depending on
#'                         on_violation setting.
#'   - bulk_fdr: The input FDR value.
#'   - null_fdr_threshold: The threshold used.
#'   - call: The matched function call.
#'
#' @examples
#' # Case 1: Bulk is non-significant (FDR = 0.15) - APPLICABLE
#' check_dds_applicability(bulk_fdr = 0.15, null_fdr_threshold = 0.05,
#'                         on_violation = "warn")
#'
#' # Case 2: Bulk is significant (FDR = 0.001) - NOT APPLICABLE
#' check_dds_applicability(bulk_fdr = 0.001, null_fdr_threshold = 0.05,
#'                         on_violation = "warn")
#'
#' @export
check_dds_applicability <- function(
    bulk_fdr,
    null_fdr_threshold = 0.05,
    on_violation = c("error", "warn", "pass")
) {

  on_violation <- match.arg(on_violation)

  # Validate inputs
  if (!is.numeric(bulk_fdr) || bulk_fdr < 0 || bulk_fdr > 1) {
    stop("bulk_fdr must be numeric in [0, 1]")
  }
  if (!is.numeric(null_fdr_threshold) || null_fdr_threshold < 0 || null_fdr_threshold > 1) {
    stop("null_fdr_threshold must be numeric in [0, 1]")
  }

  # Check applicability: bulk must be NULL (FDR > threshold)
  applicable <- bulk_fdr > null_fdr_threshold

  if (applicable) {
    # APPLICABLE case
    applicability_note <- sprintf(
      "Applicable: bulk FDR = %.4f > %.4f (null). DDS answers: 'Is this non-significant result due to dilution?'",
      bulk_fdr, null_fdr_threshold
    )
    violation_severity <- "none"
  } else {
    # NOT APPLICABLE case
    applicability_note <- sprintf(
      "Not Applicable: bulk FDR = %.4f <= %.4f (significant). DDS inapplicable; bulk effect is already detected.",
      bulk_fdr, null_fdr_threshold
    )

    if (on_violation == "error") {
      violation_severity <- "error"
      stop(applicability_note)
    } else if (on_violation == "warn") {
      violation_severity <- "warning"
      warning(applicability_note, immediate. = TRUE)
    } else {
      # on_violation == "pass"
      violation_severity <- "none"  # Suppress in pass mode
    }
  }

  return(
    list(
      applicable = applicable,
      applicability_note = applicability_note,
      violation_severity = violation_severity,
      bulk_fdr = bulk_fdr,
      null_fdr_threshold = null_fdr_threshold,
      call = match.call()
    )
  )
}
