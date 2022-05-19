#' Check for N column if not present and user wants, impute N based on user's
#' sample size. **NOTE** this will be the same value for each SNP which is not
#' necessarily correct and may cause issues down the line. N can also be
#' inputted with "ldsc", "sum", "giant" or "metal" by passing one or
#' multiple of these.
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
#' \item{\code{"ldsc"}: }{N will be computed as effective sample size:
#' Neff =(N_CAS+N_CON)*(N_CAS/(N_CAS+N_CON)) / mean((N_CAS/(N_CAS+N_CON))(N_CAS+N_CON)==max(N_CAS+N_CON))}.
#' \item{\code{"giant"}: }{N will be computed as effective sample size:
#' Neff = 2 / (1/N_CAS + 1/N_CON)}.
#' \item{\code{"metal"}: }{N will be computed as effective sample size:
#' Neff = 4 / (1/N_CAS + 1/N_CON)}.
#' }
#' @param return_list Return the \code{sumstats_dt} within a named list 
#' (default: \code{TRUE}).
#' @inheritParams compute_sample_size
#' @inheritParams read_sumstats
#' @export
#' @examples 
#' sumstats_dt <- MungeSumstats::formatted_example()
#' sumstats_dt2 <- MungeSumstats::compute_nsize(sumstats_dt=sumstats_dt,
#'                                              compute_n=10000)
compute_nsize <- function(sumstats_dt,
                          imputation_ind = FALSE,
                          compute_n = c("ldsc", "giant", "metal", "sum"),
                          standardise_headers=FALSE,
                          force_new = FALSE,
                          return_list=TRUE) {
    ## Avoid confusing BiocCheck.
    IMPUTATION_n <- IMPUTATION_Neff <- NULL;
    #### Copy data.table ####
    ## IMPORTANT!: This makes it so that any subsequent changes made to 
    ## the internal sumstats_dt variable (and the function output) will not 
    ## alter the original input data as well. 
    sumstats_dt <- data.table::copy(sumstats_dt) 
    #### Standardise the column names before continuing ####
    if (standardise_headers) {
        sumstats_dt <-
            standardise_sumstats_column_headers_crossplatform(
                sumstats_dt = sumstats_dt, 
            )[["sumstats_dt"]]
    }
    # if you want an Neff column for more than one method need to ensure
    # you can tell which is which
    append_method_name <- FALSE
    if (is.vector(compute_n)) {
        # sum makes an N column not an Neff column
        if (length(compute_n[!compute_n == "sum"]) > 1) {
              append_method_name <- TRUE
          }
    }
    for (method in compute_n) {
        compute_sample_size(
            sumstats_dt = sumstats_dt,
            method = method,
            force_new = force_new,
            append_method_name = append_method_name
        )
    }
    # if user wants information, give SNPs where Z-score calculated
    ## Evaluate the conditions under which N would have been calculate by
    ## compute_sample_size()
    if (imputation_ind &&
        is.numeric(compute_n) && compute_n != 0L &&
        (!"N" %in% names(sumstats_dt))) {
        sumstats_dt[, IMPUTATION_n := TRUE]
    }
    if (imputation_ind &&
        is.character(compute_n)) {
        sumstats_dt[, IMPUTATION_Neff := TRUE]
    }
    #### Return format ####
    if(return_list){
        return(list("sumstats_dt" = sumstats_dt))
    }else {
        return(sumstats_dt)
    } 
}
