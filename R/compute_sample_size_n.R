#' Add user supplied sample size
#'
#' @inheritParams compute_sample_size
#' @return No return
#' @keywords internal
compute_sample_size_n <- function(sumstats_dt,
                                  method,
                                  force_new = FALSE) {

    # Avoid confusing Biocheck
    N <- NULL

    if ("N" %in% names(sumstats_dt) & (!force_new)) {
        message("N already exists within sumstats_dt.")
    } else {
        if (method > 0) {
            message("Assigning N=", method, " for all SNPs.")
            sumstats_dt[, N := method]
        } else {
            msg <- "Warning: When method is an integer, must be >0."
            message(msg)
        }
    }
}
