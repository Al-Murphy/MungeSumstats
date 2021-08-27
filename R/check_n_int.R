#' Ensure that the N column is all integers
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param convert_n_int Binary, if N (the number of samples) is not an integer,
#' should this be rounded? Default is TRUE.
#' @param imputation_ind Binary Should a column be added for each imputation
#' step to show what SNPs have imputed values for differing fields. This
#' includes a field denoting SNP allele flipping (flipped). **Note**
#' these columns will be in the formatted summary statistics returned. Default
#' is FALSE.
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table :=
#' @importFrom methods is
check_n_int <- function(sumstats_dt,
                        path,
                        convert_n_int,
                        imputation_ind) {
    N <- N_tmp <- convert_n_integer <- NULL
    # Sometimes the N column is not all integers... so round it up
    col_headers <- names(sumstats_dt)
    if ("N" %in% col_headers) {
        if (!is.integer(sumstats_dt$N)) { # check if any decimal places
            msg <- paste0(
                "The sumstats N column is not all integers, this could ",
                "effect downstream analysis."
            )
            if (convert_n_int) { # if user wants to correct
                message("Ensuring that the N column is all integers.")
                msg2 <- paste0(msg, " These will be converted to integers.")
                message(msg2)
                # if user wants information, give SNPs without integer N
                if (imputation_ind) {
                      sumstats_dt[N %% 1 != 0, convert_n_integer := TRUE]
                  }
                if (is(sumstats_dt$N, "character")) {
                      suppressWarnings(sumstats_dt[, N := as.integer(N)])
                  }
                sumstats_dt[, N := round(N, 0)]
                return(list("sumstats_dt" = sumstats_dt))
            } else {
                msg2 <- paste0(msg, "These will NOT be converted to integers.")
                message(msg2)
            }
        }
    }
    return(list("sumstats_dt" = sumstats_dt))
}
