#' Vectorized Dilution Scan Across Multiple Genes
#'
#' Applies the bulknull framework to multiple genes in a single call, returning
#' a summary table with DDS, DI, and verdict for each gene. Skips bulk-significant
#' genes (FDR <= threshold) and genes with missing single-cell entries.
#'
#' @param bulk_deg data frame; differential expression results with columns
#'   `gene` (character), `beta` or `log2fc` (numeric effect size),
#'   `se` or `se_log2fc` (numeric standard error), and `fdr` or `padj`
#'   (numeric FDR). Additional columns are ignored.
#' @param d_sc_values named numeric vector; standardized effect sizes keyed by gene (primary input)
#' @param sc_zscores named numeric vector; z-scores keyed by gene (backward compatibility).
#'   Genes not in either vector are skipped with a note.
#' @param cluster_fraction numeric; fraction of total cells in target cluster
#'   (constant across all genes).
#' @param condition_fraction numeric; fraction of cells in case condition
#'   (constant across all genes).
#' @param r numeric; relative expression parameter. Default 1.
#' @param n_cluster numeric; number of cells in target cluster
#'   (constant across all genes). Required if sc_zscores supplied.
#' @param n_eff_basis character; "cell" or "donor" basis for z-score conversion.
#' @param n_donors numeric; number of donors (required if n_eff_basis = "donor").
#' @param null_fdr_threshold numeric; FDR threshold for applicability gate
#'   (default 0.05). Genes with FDR <= threshold are marked "FAILED".
#' @param alpha numeric; significance level for DI (default 0.05).
#' @param sided character; "one" or "two" for one-sided or two-sided test
#'   (default "one").
#'
#' @return A data frame with columns:
#'   \item{gene}{Gene name from bulk_deg}
#'   \item{bulk_beta}{Effect size from bulk_deg}
#'   \item{bulk_se}{Standard error from bulk_deg}
#'   \item{bulk_fdr}{FDR from bulk_deg}
#'   \item{sc_zscore}{Single-cell z-score (NA if not found)}
#'   \item{mu_dilution}{Expected bulk z-score under dilution hypothesis}
#'   \item{dds}{Dilution Discordance Score}
#'   \item{di}{Diagnosability Index}
#'   \item{dds_interpretation}{Categorical interpretation of DDS}
#'   \item{di_interpretation}{Categorical interpretation of DI}
#'   \item{verdict}{Final interpretive verdict}
#'   \item{gate_status}{PASSED or FAILED (applicability gate status)}
#'
#' @details
#' The function iterates over rows of bulk_deg, skipping:
#' - Genes with bulk FDR <= null_fdr_threshold (gate fails)
#' - Genes without an entry in sc_zscores
#' - Genes with zero variance in the single-cell data (sc_zscore = 0)
#'
#' For applicable genes, it calls bulknull() internally and extracts key results.
#' Returns all processed genes in output, with NA or "SKIPPED" for missing data.
#'
#' @examples
#' # Example: scan a small DEG table with mock z-scores
#' bulk_deg <- data.frame(
#'   gene = c("GENE1", "GENE2", "GENE3"),
#'   beta = c(-0.084, 0.020, 0.150),
#'   se = c(0.063, 0.050, 0.100),
#'   fdr = c(0.450, 0.05, 0.001)
#' )
#' sc_zscores <- c(GENE1 = 2.077, GENE2 = 1.5, GENE3 = 0.5)
#'
#' result <- dilution_scan(
#'   bulk_deg = bulk_deg,
#'   sc_zscores = sc_zscores,
#'   cluster_fraction = 0.0695,
#'   condition_fraction = 0.631,
#'   n_cluster = 1332
#' )
#' print(result)
#'
#' @export
dilution_scan <- function(bulk_deg, d_sc_values = NULL, sc_zscores = NULL,
                          cluster_fraction, condition_fraction, r = 1, n_cluster = NULL,
                          n_eff_basis = c("cell", "donor"), n_donors = NULL,
                          null_fdr_threshold = 0.05, alpha = 0.05,
                          sided = c("one", "two")) {

  sided <- match.arg(sided)
  n_eff_basis <- match.arg(n_eff_basis)

  # Validate input
  if (is.null(d_sc_values) && is.null(sc_zscores)) {
    stop("Must supply either d_sc_values or sc_zscores")
  }

  # Identify column names (flexible naming for FDR, effect size, SE)
  fdr_col <- if ("fdr" %in% names(bulk_deg)) "fdr" else if ("padj" %in% names(bulk_deg)) "padj" else stop("bulk_deg must have 'fdr' or 'padj' column")
  beta_col <- if ("beta" %in% names(bulk_deg)) "beta" else if ("log2fc" %in% names(bulk_deg)) "log2fc" else stop("bulk_deg must have 'beta' or 'log2fc' column")
  se_col <- if ("se" %in% names(bulk_deg)) "se" else if ("se_log2fc" %in% names(bulk_deg)) "se_log2fc" else stop("bulk_deg must have 'se' or 'se_log2fc' column")
  gene_col <- if ("gene" %in% names(bulk_deg)) "gene" else if ("rownames" %in% names(bulk_deg)) "rownames" else stop("bulk_deg must have 'gene' column")

  # Initialize output
  results <- list()

  for (i in seq_len(nrow(bulk_deg))) {
    gene <- bulk_deg[[i, gene_col]]
    beta <- bulk_deg[[i, beta_col]]
    se <- bulk_deg[[i, se_col]]
    fdr <- bulk_deg[[i, fdr_col]]

    # Get d_sc or sc_zscore (d_sc takes priority)
    d_sc <- NULL
    sc_zscore <- NULL

    if (!is.null(d_sc_values) && gene %in% names(d_sc_values) && !is.na(d_sc_values[[gene]])) {
      d_sc <- d_sc_values[[gene]]
    } else if (!is.null(sc_zscores) && gene %in% names(sc_zscores) && !is.na(sc_zscores[[gene]])) {
      sc_zscore <- sc_zscores[[gene]]
    } else {
      results[[i]] <- list(
        gene = gene,
        bulk_beta = beta,
        bulk_se = se,
        bulk_fdr = fdr,
        sc_zscore = NA,
        d_sc = NA,
        mu_dilution = NA,
        dds = NA,
        di = NA,
        dds_interpretation = "SKIPPED (missing d_sc/z-score)",
        di_interpretation = "SKIPPED (missing d_sc/z-score)",
        verdict = "SKIPPED: Single-cell values not provided",
        gate_status = "SKIPPED"
      )
      next
    }

    # Check for zero signal
    if ((!is.null(sc_zscore) && sc_zscore == 0) || (!is.null(d_sc) && d_sc == 0)) {
      results[[i]] <- list(
        gene = gene,
        bulk_beta = beta,
        bulk_se = se,
        bulk_fdr = fdr,
        sc_zscore = sc_zscore,
        d_sc = d_sc,
        mu_dilution = NA,
        dds = NA,
        di = NA,
        dds_interpretation = "SKIPPED (zero signal)",
        di_interpretation = "SKIPPED (zero signal)",
        verdict = "SKIPPED: Zero single-cell signal",
        gate_status = "SKIPPED"
      )
      next
    }

    # Apply bulknull gate and compute diagnostics
    tryCatch({
      bulknull_result <- bulknull(
        bulk_beta = beta,
        bulk_se = se,
        bulk_fdr = fdr,
        d_sc = d_sc,
        sc_zscore = sc_zscore,
        cluster_fraction = cluster_fraction,
        condition_fraction = condition_fraction,
        r = r,
        n_cluster = n_cluster,
        n_eff_basis = n_eff_basis,
        n_donors = n_donors,
        null_fdr_threshold = null_fdr_threshold,
        on_violation = "pass",
        alpha = alpha,
        sided = sided
      )

      results[[i]] <- list(
        gene = gene,
        bulk_beta = beta,
        bulk_se = se,
        bulk_fdr = fdr,
        sc_zscore = sc_zscore,
        d_sc = bulknull_result$d_sc,
        mu_dilution = bulknull_result$mu_dilution,
        dds = bulknull_result$.dds_full$dds_score,
        di = bulknull_result$di,
        dds_interpretation = bulknull_result$dds_interpretation,
        di_interpretation = bulknull_result$di_interpretation,
        verdict = bulknull_result$verdict,
        gate_status = bulknull_result$gate_status
      )
    }, error = function(e) {
      results[[i]] <<- list(
        gene = gene,
        bulk_beta = beta,
        bulk_se = se,
        bulk_fdr = fdr,
        sc_zscore = sc_zscore,
        d_sc = d_sc,
        mu_dilution = NA,
        dds = NA,
        di = NA,
        dds_interpretation = "ERROR",
        di_interpretation = "ERROR",
        verdict = paste("ERROR:", e$message),
        gate_status = "ERROR"
      )
    })
  }

  # Convert to data frame (build by column to handle NULL values correctly)
  if (length(results) == 0) {
    return(data.frame())
  }

  # Helper function to extract value, handling NULL -> NA
  extract_with_na <- function(lst, key) {
    vals <- lapply(lst, function(x) {
      val <- x[[key]]
      if (is.null(val)) NA else val
    })
    # Try to coerce to numeric, fall back to character if not numeric
    tryCatch(
      as.numeric(unlist(vals)),
      warning = function(w) unlist(vals)
    )
  }

  # Extract columns from all results
  data.frame(
    gene = sapply(results, "[[", "gene"),
    bulk_beta = extract_with_na(results, "bulk_beta"),
    bulk_se = extract_with_na(results, "bulk_se"),
    bulk_fdr = extract_with_na(results, "bulk_fdr"),
    sc_zscore = extract_with_na(results, "sc_zscore"),
    d_sc = extract_with_na(results, "d_sc"),
    mu_dilution = extract_with_na(results, "mu_dilution"),
    dds = extract_with_na(results, "dds"),
    di = extract_with_na(results, "di"),
    dds_interpretation = sapply(results, "[[", "dds_interpretation"),
    di_interpretation = sapply(results, "[[", "di_interpretation"),
    verdict = sapply(results, "[[", "verdict"),
    gate_status = sapply(results, "[[", "gate_status"),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
