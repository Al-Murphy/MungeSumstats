#' Compute (effective) sample size
#'
#' Computes sample sum (as new column "N") or
#' effective sample size (ESS) (as new column "Neff").
#' Computing ESS is important as it takes into account
#' the proportion of cases to controls (i.e. class imbalance) so as not to
#' overestimate your statistical power.
#'
#' There are many different formulas for calculating ESS,
#' but LDSC is probably the best method available here, as it
#'  doesn't assume that the proportion of controls:cases
#'  is 2:1 (as in GIANT) or 4:1 (as in METAL).
#'
#' @param sumstats_dt Summary statistics data.table.
#' @param method Method for computing (effective) sample size.
#'
#' \itemize{
#'
#' \item{"ldsc" : \cr}{
#' \eqn{Neff = (N_CAS+N_CON) * (N_CAS/(N_CAS+N_CON)) /
#'  mean((N_CAS/(N_CAS+N_CON))[(N_CAS+N_CON)==max(N_CAS+N_CON)]))}\cr
#' \href{https://github.com/bulik/ldsc/issues/95}{
#' bulik/ldsc GitHub Issue}
#' \href{https://github.com/bulik/ldsc/blob/aa33296abac9569a6422ee6ba7eb4b902422cc74/munge_sumstats.py#L321}{
#' bulik/ldsc GitHub code}
#' }
#'
#' \item{"giant" : \cr}{
#' \eqn{Neff = 2 / (1/N_CAS + 1/N_CON)}\cr
#' \href{https://www.nature.com/articles/nprot.2014.071}{
#' Winkler et al. 2014, Nature}
#' }
#'
#' \item{"metal" : \cr}{
#' \eqn{Neff = 4 / (1/N_CAS + 1/N_CON)}\cr
#' \href{https://pubmed.ncbi.nlm.nih.gov/20616382/}{
#' Willer et al. 2010, Bioinformatics}
#' }
#'
#' \item{"sum" : \cr}{
#' \eqn{N = N_CAS + N_CON}\cr
#' Simple summation of cases and controls
#'  that does not account for class imbalance.
#'  }
#'
#' \item{"\\<integer\\>" : \cr}{
#' \code{N = \\<integer\\>}\cr
#' If method is a positive integer, it will be used as N
#' for every row.
#' }
#'
#' }
#'
#' @param force_new If "Neff" (or "N") already exists in \code{sumstats_dt},
#' replace it with the recomputed version.
#' @param append_method_name should Neff column have an indicator to explain the 
#' method that makes it., Default is FALSE unless multiple methods are passed
#'
#' @return A data.table with a new column "Neff" or "N"
#' @keywords internal
compute_sample_size <- function(sumstats_dt,
                                method = c("ldsc", "giant", "metal", "sum"),
                                force_new = FALSE,
                                append_method_name=FALSE) {
    if (is.character(method)) {
        #### Compute Neff/N ###
        compute_sample_size_neff(
            sumstats_dt = sumstats_dt,
            method = method,
            force_new = force_new,
            append_method_name=append_method_name
        )
    } else if (is.numeric(method)) {
        #### Add user-supplied N ###
        compute_sample_size_n(
            sumstats_dt = sumstats_dt,
            method = method,
            force_new = force_new
        )
    } else {
        all_method <- c("ldsc", "giant", "metal", "sum", "<integer>")
        msg <- paste0(
            "Warning: method must be one of:\n",
            paste0(" - ", all_method, collapse = "\n")
        )
        message(msg)
    }
}
