#' Check numeric columns
#' 
#' Checks for any columns that should be numeric, 
#' and ensures that they are indeed numeric.
#' @param sumstats_dt Summary stats with column names already standardised by
#'  \link[MungeSumstats]{format_sumstats}.
#' @param cols Names of columns that should be numeric. 
#' If any of these columns are not actually present in \code{sumstats_dt},
#' they will be skipped.
#' @returns sumstats_dt 
#' 
#' @keywords internal
#' @importFrom data.table :=
check_numeric <- function(sumstats_dt,
                          cols = c("P","SE","FRQ","MAF","BETA")){
    cols <- cols[cols %in% names(sumstats_dt)]
    if(length(cols)>0){
        sumstats_dt[,(cols):=lapply(.SD, as.numeric), .SDcols=cols]   
    }
    return(sumstats_dt)
}