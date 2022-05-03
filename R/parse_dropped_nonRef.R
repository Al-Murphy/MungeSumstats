#' Parse number of SNPs dropped due to being in the ref genome
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_dropped_nonRef <- function(l){
    line <- grep("are not on the reference genome",l,value = TRUE)[1]
    if(is.na(line)) return(NA)
    as.integer(trimws(gsub(",","",strsplit(line," SNPs")[[1]][1])))
}