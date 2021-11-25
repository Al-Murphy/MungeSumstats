#' Parse number of SNPs with p-values >1
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_pval_large <- function(l){
    line <- grep(
        "p-values are >1 which LDSC/MAGMA may not be able to handle",l,
        value = TRUE)[1]
    as.integer(trimws(gsub(",","",strsplit(line," ")[[1]][1])))
}