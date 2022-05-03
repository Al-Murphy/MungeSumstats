#' Parse number of SNPs with p-values <0
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_pval_neg <- function(l){
    line <- grep(
        "p-values are <0 which LDSC/MAGMA may not be able to handle",l,
        value = TRUE)[1]
    if(is.na(line)) return(NA)
    as.integer(trimws(gsub(",","",strsplit(line," ")[[1]][1])))
}