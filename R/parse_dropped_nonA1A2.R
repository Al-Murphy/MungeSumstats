#' Parse number of SNPs dropped due to not matching the ref genome A1 or A2
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_dropped_nonA1A2 <- function(l){
    line <- grep("neither A1 nor A2 match the reference genome",l,
                 value = TRUE)[1]
    as.integer(trimws(gsub("There are|,","",strsplit(line," SNPs")[[1]][1])))
}