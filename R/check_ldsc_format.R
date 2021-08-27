#' Ensures that parameters are compatible with LDSC format
#'
#' Format summary statistics for direct input to
#' Linkage Disequilibrium SCore (LDSC) regression without the need
#' to use their \code{munge_sumstats.py} script first.
#'
#' \href{https://github.com/bulik/ldsc/wiki/Summary-Statistics-File-Format}{LDSC documentation}.
#'
#' @inheritParams format_sumstats
#' @inheritParams check_zscore
#'
#' @return Formatted summary statistics
#' @source \href{https://github.com/bulik/ldsc}{LDSC GitHub}
check_ldsc_format <- function(sumstats_dt, ldsc_format, convert_n_int,
                              allele_flip_check, compute_z, compute_n) {
    # LDSC Requires SNP,N,Z,A1,A2
    # SNP,A1,A2 is enforced by MungeSumstats automatically
    if (ldsc_format) {
        message("Ensuring parameters comply with LDSC format.")
        # check if columns already exist
        z_present <- "Z" %in% names(sumstats_dt)
        n_present <- "N" %in% names(sumstats_dt)
        if (!convert_n_int) {
            message("Setting `convert_n_int=TRUE` to comply with LDSC format.")
            convert_n_int <- TRUE
        }
        if (!allele_flip_check) {
            message("Setting `allele_flip_check=TRUE` to comply with LDSC format.")
            allele_flip_check <- TRUE
        }
        if (!compute_z && !z_present) {
            message("Setting `compute_z=TRUE` to comply with LDSC format.")
            compute_z <- TRUE
        }
        n_msg <- paste0(
            "LDSC requires an N column but your dataset doesn't ",
            "appear to have one. You can impute an N value for ",
            "all your SNPs\nby setting `compute_n` to this value but",
            " note this is may not be correct and may lead to ",
            "different results from LDSC\nthan if the true N per SNP",
            " was known."
        )
        if (compute_n == 0L && !n_present) {
              stop(n_msg)
          }
    }
    return(list(
        ldsc_format = ldsc_format,
        convert_n_int = convert_n_int,
        allele_flip_check = allele_flip_check,
        compute_z = compute_z
    ))
}
