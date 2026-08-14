#' DPP9/M8 Example Dataset
#'
#' A list containing measurements for DPP9 (diphosphonucleotidase 9) expression
#' in the M8 (inflammatory macrophages) subcluster from IPF lung tissue.
#'
#' This is a real data example from the bulknull application: a bulk-level null
#' result (non-significant effect) paired with strong single-cell signal in a
#' small subcluster. Used throughout vignettes and examples.
#'
#' @format A list with 12 elements:
#'   \describe{
#'     \item{bulk_beta}{Numeric. Bulk effect size (-0.0297)}
#'     \item{bulk_se}{Numeric. Standard error of bulk effect (0.0595)}
#'     \item{bulk_fdr}{Numeric. FDR-adjusted p-value from meta-analysis (0.821)}
#'     \item{sc_zscore}{Numeric. Single-cell module z-score (2.077)}
#'     \item{cluster_fraction}{Numeric. Fraction of tissue in M8 (0.0695)}
#'     \item{condition_fraction}{Numeric. Fraction of M8 from IPF (0.631)}
#'     \item{n_cluster}{Integer. Number of cells in M8 (1332)}
#'     \item{gene}{Character. Gene name ("DPP9")}
#'     \item{cluster}{Character. Cluster identifier ("M8")}
#'     \item{cluster_name}{Character. Cluster description ("Inflammatory Macrophages")}
#'     \item{cohort}{Character. Disease cohort ("IPF")}
#'     \item{description}{Character. Dataset description}
#'   }
#'
#' @details
#' The DPP9/M8 case exemplifies the framework's utility: bulk RNA-seq from
#' IPF lung tissue shows no significant DPP9 effect (FDR = 0.821), but
#' single-cell analysis reveals strong inflammasome-module enrichment in M8.
#' The DDS computes whether the bulk null reflects true absence or dilution
#' of this M8 signal by the dominant other cell types.
#'
#' **Result:** DDS ~ 0.48 (ambiguous), suggesting both null and dilution
#' remain plausible. The Diagnosability Index (DI ~ 0.065) reveals that the
#' study design is fundamentally underpowered for dilution detection due to
#' the small M8 fraction (6.95% of lung), independent of effect size.
#'
#' @usage data(dpp9_m8)
#'
#' @examples
#' data(dpp9_m8)
#' str(dpp9_m8)
#'
#' # Compute Dilution Discordance Score
#' result <- dilution_score(
#'   bulk_beta = dpp9_m8$bulk_beta,
#'   bulk_se = dpp9_m8$bulk_se,
#'   bulk_fdr = dpp9_m8$bulk_fdr,
#'   sc_zscore = dpp9_m8$sc_zscore,
#'   cluster_fraction = dpp9_m8$cluster_fraction,
#'   condition_fraction = dpp9_m8$condition_fraction,
#'   n_cluster = dpp9_m8$n_cluster
#' )
#' print(result$dds_score)
"dpp9_m8"
