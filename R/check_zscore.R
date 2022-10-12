#' Check for Z-score column
#'
#' The following ensures that a Z-score column is present.
#' The Z-score formula we used here is a R implementation of the formula
#' used in \href{https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/munge_sumstats.py#L363}{LDSC's munge_sumstats.py}:
#'
#' \code{np.sqrt(chi2.isf(P, 1))}
#'
#' The R implementation is adapted from the \code{GenomicSEM::munge} function,
#' after optimizing for speed using \code{data.table}:
#'
#' \code{sumstats_dt[,Z:=sign(BETA)*sqrt(stats::qchisq(P,1,lower=FALSE))]}
#'
#' \emph{NOTE}: \code{compute_z} is set to \code{TRUE} by 
#' default to ensure standardisation
#' of the "Z" column (which can be computed differently in different datasets).
#'
#' @return \code{list("sumstats_dt"=sumstats_dt)}
#'
#' @param sumstats_dt data table obj of the summary statistics file for the
#' GWAS.
#' @param imputation_ind Binary Should a column be added for each imputation
#' step to show what SNPs have imputed values for differing fields. This
#' includes a field denoting SNP allele flipping (flipped). **Note**
#' these columns will be in the formatted summary statistics returned. Default
#' is FALSE.
#' @param compute_z Whether to compute Z-score column. Default is FALSE. This 
#' can be computed from Beta and SE with (Beta/SE) or P 
#' (Z:=sign(BETA)*sqrt(stats::qchisq(P,1,lower=FALSE))).
#' **Note** that imputing the Z-score from P for every SNP will not be
#' perfectly correct and may result in a loss of power. This should only be done
#' as a last resort. Use 'BETA' to impute by BETA/SE and 'P' to impute by SNP
#' p-value.
#' @param force_new_z When a "Z" column already exists, it will be used by
#' default. To override and compute a new Z-score column from P set
#' \code{force_new_z=TRUE}.
#' @param standardise_headers Run
#' \code{standardise_sumstats_column_headers_crossplatform} first.
#' @param mapping_file MungeSumstats has a pre-defined column-name mapping file
#' which should cover the most common column headers and their interpretations.
#' However, if a column header that is in youf file is missing of the mapping we
#' give is incorrect you can supply your own mapping file. Must be a 2 column
#' dataframe with column names "Uncorrected" and "Corrected". See
#' data(sumstatsColHeaders) for default mapping and necessary format.
#'
#' @keywords internal
#' @importFrom stats qchisq
check_zscore <- function(sumstats_dt,
                         imputation_ind,
                         compute_z = 'BETA',
                         force_new_z = FALSE,
                         standardise_headers = FALSE,
                         mapping_file) {
    ## Set variables to be used in in place data.table functions to NULL
    ## to avoid confusing BiocCheck.
    Z <- BETA <- SE <- P <- IMPUTATION_z_score_beta_se <- 
      IMPUTATION_z_score_p <- NULL

    if (standardise_headers) {
        sumstats_dt <-
            standardise_sumstats_column_headers_crossplatform(
                sumstats_dt = sumstats_dt,
                mapping_file = mapping_file
            )[["sumstats_dt"]]
    }

    if (!isFALSE(compute_z)) {
        # message("Checking Z-score.")
        col_headers <- names(sumstats_dt)
        if ("Z" %in% col_headers && (!force_new_z)) {
            message("Keeping existing Z-score column.")
        } else if(toupper(compute_z)=='BETA'){
            if ("BETA" %in% col_headers && "SE" %in% col_headers){
              message(
                "Computing Z-score from BETA ans SE using formula:",
                " `BETA/SE`"
              )
              # ensure BETA, P are numeric
              sumstats_dt[, BETA := as.numeric(BETA)]
              sumstats_dt[, SE := as.numeric(SE)]
              sumstats_dt[, Z := BETA/SE]
              # if user wants information, give SNPs where Z-score calculated
              if (imputation_ind) {
                sumstats_dt[, IMPUTATION_z_score_beta_se := TRUE]
              }
            }
            else{
              msg <- paste0("Can't compute Z-score from BETA and SE as both ",
                            "aren't in the sumstats.\nPlease choose another ",
                            "option for `compute_z`")
              stop(msg)
            }
          
        } else{ #P-value
            message(
                "Computing Z-score from P using formula:",
                " `sign(BETA)*sqrt(stats::qchisq(P,1,lower=FALSE)`"
            )
            # ensure BETA, P are numeric
            sumstats_dt[, BETA := as.numeric(BETA)]
            sumstats_dt[, P := as.numeric(P)]
            sumstats_dt[, Z := sign(BETA) *
                sqrt(stats::qchisq(P, 1, lower = FALSE))]
            # if user wants information, give SNPs where Z-score calculated
            if (imputation_ind) {
                sumstats_dt[, IMPUTATION_z_score_p := TRUE]
            }
        }
    }
    return(list("sumstats_dt" = sumstats_dt))
}
