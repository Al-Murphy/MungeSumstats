#' Parse number of SNPs flipped to align with the ref genome
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_flipped <- function(l){
    line <- grep("where A1 doesn't match the reference genome",l,
                 value = TRUE)[1]
    if(is.na(line)) return(NA)
    as.integer(trimws(gsub("There are|,","",strsplit(line," SNPs")[[1]][1])))
}