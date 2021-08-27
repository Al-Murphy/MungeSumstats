#' Check for N column if not present and user wants, impute N based on user's
#' sample size. **NOTE** this will be the same value for each SNP which is not
#' necessarily correct and may cause issues down the line.
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
#' @param compute_n How to compute per-SNP sample size (new column "N").
#' \itemize{
#' \item{\code{0}: }{N will not be computed.}
#' \item{\code{>0}: }{If any number >0 is provided,
#' that value will be set as N for every row.
#' **Note**: Computing N this way is incorrect and should be avoided
#' if at all possible.}
#' \item{\code{"sum"}: }{N will be computed as:
#' cases (N_CAS) + controls (N_CON), so long as both columns are present}.
#' \item{\code{"effective"}: }{N will be computed as effective sample size:
#' cases (N_CAS) + controls (N_CON), so long as both columns are present}.
#' } 
#'
#' @keywords internal
compute_nsize <- function(sumstats_dt,
                          imputation_ind = FALSE,
                          compute_n = c("ldsc", "sum"),
                          force_new = FALSE) {
    ## Avoid confusing BiocCheck.
    IMPUTATION_n <- NULL

    for (method in compute_n) {
        compute_sample_size(
            sumstats_dt = sumstats_dt,
            method = method,
            force_new = force_new
        )
    }
    # if user wants information, give SNPs where Z-score calculated
    ## Evaluate the conditions under which N would have been calculate by
    ## compute_sample_size()
    if (imputation_ind &
        is.numeric(compute_n) &
        (!"N" %in% names(sumstats_dt))) {
        sumstats_dt[, IMPUTATION_n := TRUE]
    }
    return(list("sumstats_dt" = sumstats_dt))
}
