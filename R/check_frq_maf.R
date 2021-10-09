#' Check that FRQ column refers to minor/effect allele frequency not major
#'
#' @inheritParams format_sumstats
#' @return sumstats_dt, the modified summary statistics data table object
#' @keywords internal
check_frq_maf <- function(sumstats_dt, frq_is_maf) {
    ## Set variables to be used in in place data.table functions to NULL
    ## to avoid confusing BiocCheck.
    FRQ <- .N <- NULL
    col_headers <- names(sumstats_dt)
    if ("FRQ" %in% col_headers) {
        # get proportion of SNPs with FRQ>0.5, this will mean major for
        # bi-allelic SNPs but not major allele frq for
        # non bi-allelic SNPs may be lower
        num_major <- sumstats_dt[FRQ > 0.5, .N, ]
        # only continue if there are some
        if (num_major > 0) {
            per_major <- round(num_major / nrow(sumstats_dt) * 100, 1)
            # get mappings for message
            frq_choices <-
                paste0(
                    sumstatsColHeaders[
                        sumstatsColHeaders$Corrected == "FRQ", ]$Uncorrected, 
                    collapse = ", ")
            msg <-
                paste0(
                    formatC(num_major, big.mark = ",", format = "fg"),
                    " SNPs (", per_major, "%) have FRQ values > 0.5. ",
                    "Conventionally the FRQ column is intended to show the",
                    " minor/effect allele frequency.\nThe FRQ column was",
                    " mapped from one of the following from the inputted ",
                    " summary statistics file:\n", frq_choices
                )
            message(msg)
            # if frq is minor allele frequency is set to TRUE 
            # just accept that it is and don't rename FRQ. 
            # If FALSE and  there are some SNPS with FRQ>0.5 then do
            if (isFALSE(frq_is_maf)) {
                msg <- paste0(
                    "As frq_is_maf=FALSE, the FRQ column will be ",
                    "renamed MAJOR_ALLELE_FRQ to differentiate the",
                    " values from \nminor/effect allele frequency."
                )
                message(msg)
                setnames(sumstats_dt, "FRQ", "MAJOR_ALLELE_FRQ")
            } else { # frq_is_maf =TRUE, i.e. don't rename
                msg <- paste0(
                    "As frq_is_maf=TRUE, the FRQ column will not be ",
                    "renamed. If the FRQ values were intended to ",
                    "represent major allele frequency,\nset ",
                    "frq_is_maf=FALSE to rename the column as ",
                    "MAJOR_ALLELE_FRQ and differentiate it",
                    " from minor/effect allele frequency."
                )
                message(msg)
            }
        }
    }
    return(sumstats_dt)
}
