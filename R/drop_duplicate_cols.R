#' Drop duplicate columns
#' 
#' Drop columns with identical names (if any exist) within a data.table.
#' @param dt data.table
#' @keywords internal
#' @return Null output
drop_duplicate_cols <- function(dt){
    dups <- which(duplicated(names(dt)))
    if(length(dups)>0){
        messager("Dropping",length(dups),"duplicate columns.")
        dt[,  which(duplicated(names(dt))):= NULL] 
    } 
}
