#' Ensure A1 and A2 are upper case
#'
#' @inheritParams format_sumstats
#' @param log_files list of log file locations
#' @return list containing sumstats_dt, the modified summary statistics data
#' table object and the log file list
#' @keywords internal
make_allele_upper <- function(sumstats_dt, log_files) {
    A1 <- A2 <- NULL
    col_headers <- names(sumstats_dt)
    
    if ("A1" %in% col_headers) {
        message("Checking A1 is uppercase")
        sumstats_dt[, A1:=toupper(A1)]
    }
    if ("A2" %in% col_headers) {
        message("Checking A2 is uppercase")
        sumstats_dt[, A2:=toupper(A2)]
    }
    
    return(list("sumstats_dt" = sumstats_dt, "log_files" = log_files))
}
