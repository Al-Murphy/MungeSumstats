#' Parse number of SNPs dropped due to being below the INFO threshold
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_dropped_INFO <- function(l){
    line <- grep("SNPs are below the INFO threshold of",l,
                 value = TRUE)[1]
    if(is.na(line)) return(NA)
    as.integer(trimws(gsub("There are|,","",strsplit(line," SNPs")[[1]][1])))
}