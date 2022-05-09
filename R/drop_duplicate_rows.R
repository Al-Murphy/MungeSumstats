#' Drop duplicate rows
#' 
#' Drop rows with duplicate values across all columns. 
#' @param dt data.table
#' @param verbose Print messages.
#' @keywords internal
#' @return Filtered \code{dt}. 
drop_duplicate_rows <- function(dt,
                                verbose=TRUE){
    nrows_old <- nrow(dt)
    dt <- unique(dt)
    nrows_new <- nrow(dt)
    dropped <-  nrows_old-nrows_new
    if(dropped>0){
        messager("Dropped",formatC(dropped,big.mark = ","),
                 "duplicate rows.",v=verbose)
    }
    return(dt)
}
