#' Remove empty columns
#' 
#' Remote columns that are empty or contain all the same values in a data.table.
#' @inheritParams format_sumstats
#' @inheritParams check_empty_cols
#' @keywords internal
#' @returns Null output.
remove_empty_cols <- function(sumstats_dt, 
                              sampled_rows=NULL,
                              verbose=TRUE){
    #### Check for empty columns ####
    messager("Checking for empty columns.",v=verbose)
    empty_cols <- check_empty_cols(
        sumstats_dt = sumstats_dt,
        sampled_rows = sampled_rows,
        verbose = FALSE
    )
    #### Remove empty columns #####
    if(length(empty_cols)>0) {
        messager("Removing",length(empty_cols),"empty columns.",v=verbose)
        sumstats_dt[,(names(empty_cols)):=NULL] 
    }
}
