#' Parse number of SNPs dropped due to being duplicates
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_dropped_duplicates <- function(l){
    line <- grep("are duplicated in the sumstats file",l,
                 value = TRUE)[1]
    if(is.na(line)) return(NA)
    as.integer(trimws(gsub(",","",strsplit(line," ")[[1]][1])))
}