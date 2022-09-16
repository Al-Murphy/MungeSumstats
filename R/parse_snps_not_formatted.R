#' Parse number of SNPs not correctly formatted
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_snps_not_formatted <- function(l){
    line <- grep("SNP IDs are not correctly formatted.",l,
                 value = TRUE)[1]
    if(is.na(line)) return(NA)
    as.integer(trimws(gsub("There are|,","",strsplit(line," SNPs")[[1]][1])))
}