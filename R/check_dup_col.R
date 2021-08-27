#' Ensure that no columns are duplicated
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
check_dup_col <- function(sumstats_dt, path) {
    message("Checking for duplicate columns.")
    col_headers <- names(sumstats_dt)
    if (sum(duplicated(col_headers)) > 0) {
        msg <- paste0(
            "There are ", sum(duplicated(col_headers)),
            " duplicated columns which will be removed."
        )
        message(msg)
        notDup <- which(!duplicated(col_headers))
        sumstats_dt <- sumstats_dt[, notDup, with = FALSE]
        return(list("sumstats_dt" = sumstats_dt))
    } else {
        return(list("sumstats_dt" = sumstats_dt))
    }
}
