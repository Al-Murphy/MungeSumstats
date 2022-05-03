#' Unlist a data.table
#' 
#' Identify columns that are lists and turn them into vectors.
#' @param dt data.table
#' @keywords internal
#' @importFrom data.table .SD
#' @importFrom methods is
unlist_dt <- function(dt) {
    .SD <- NULL
    cols <- names(dt)[sapply(dt, methods::is,"list")]
    if(length(cols)>0){
        messager("Unlisting",length(cols),"columns.")
        dt[ , (cols) := lapply(.SD,unlist), .SDcols = cols]
    } 
}
