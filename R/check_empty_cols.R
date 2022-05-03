#' Check for empty columns
#'
#' Empty columns contain only ".", NA, or 0
#' @inheritParams format_sumstats
#' @param sampled_rows First N rows to sample.
#' Set \code{NULL} to use full \code{sumstats_file}.
#' when determining whether cols are empty.
#' 
#' @param verbose Print messages.
#' 
#' @return empty_cols
#' @keywords internal
#' @importFrom utils head
check_empty_cols <- function(sumstats_dt,
                             sampled_rows = NULL, 
                             verbose = TRUE) {
    if (is.null(sampled_rows)) {
        sampled_rows <- nrow(sumstats_dt)
    } else {
        sampled_rows <- min(sampled_rows, nrow(sumstats_dt))
    }
    empty_cols <- vapply(
        colnames(sumstats_dt), function(x) {
            dt <- utils::head(sumstats_dt, sampled_rows)
        (sum(unlist(dt[[x]]) != ".") == 0) |
            (sum(!is.na(unlist(dt[[x]]))) == 0) |
            (sum(unlist(dt[[x]]) != 0) == 0)
    },
    FUN.VALUE = logical(1)
    )
    empty_cols <- empty_cols[empty_cols]
    messager(length(empty_cols), "empty column(s) detected.",
             v=verbose)
    return(empty_cols)
}
