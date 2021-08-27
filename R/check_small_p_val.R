#' Ensure that the p values are not 5e-324 or lower, if so set to 0
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param convert_small_p Binary, should p-values < 5e-324 be converted to 0?
#' Small p-values pass the R limit and can cause errors with LDSC/MAGMA and
#' should be converted. Default is TRUE.
#' @param imputation_ind Binary Should a column be added for each imputation
#' step to show what SNPs have imputed values for differing fields. This
#' includes a field denoting SNP allele flipping (flipped). **Note**
#' these columns will be in the formatted summary statistics returned. Default
#' is FALSE.
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table :=
check_small_p_val <- function(sumstats_dt, path, convert_small_p,
                              imputation_ind) {
    P <- convert_small_p_0 <- NULL
    # Sometimes the N column is not all integers... so round it up
    col_headers <- names(sumstats_dt)
    if ("P" %in% col_headers) {
        # get smallest p-val - seems to change to character if < xe-300
        char_check <- FALSE
        num_check <- FALSE
        if (is.numeric(sumstats_dt$P)) {
            if (min(sumstats_dt$P, na.rm = TRUE) <= 5e-324) {
                  num_check <- TRUE
              }
        } else { # char check
            max_minus_power <- max(as.numeric(gsub(".*-", "", sumstats_dt$P)))
            if (max_minus_power >= 324) {
                  char_check <- TRUE
              }
        }
        if (char_check | num_check) { # check if any smaller or equal to 5e-324 limit
            msg <- paste0(
                "There are existing p-values as low as 5e-324 which ",
                "LDSC/MAGMA may not be able to handle. "
            )
            if (convert_small_p) { # if user wants to correct
                msg2 <- paste0(msg, "These will be converted to 0.")
                message(msg2)
                sumstats_dt[, P := as.numeric(as.character(P))]
                # if users want edited snps, return information
                if (imputation_ind) {
                      sumstats_dt[P == 0, convert_small_p_0 := TRUE]
                  }

                return(list("sumstats_dt" = sumstats_dt))
            } else {
                msg2 <- paste0(msg, "These will NOT be converted.")
                message(msg2)
            }
        }
    }
    return(list("sumstats_dt" = sumstats_dt))
}
