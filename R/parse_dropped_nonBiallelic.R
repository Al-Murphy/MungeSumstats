#' Parse number of SNPs dropped due to not being bi-allelic
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_dropped_nonBiallelic <- function(l){
    line <- grep("SNPs are non-biallelic",l,
                 value = TRUE)[1]
    as.integer(trimws(gsub(",","",strsplit(line," ")[[1]][1])))
}