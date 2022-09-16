#' Parse number/percent of SNPs with FREQ values >0.5
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_snps_freq_05 <- function(l,
                               percent=FALSE){
    line <- grep("have FRQ values > 0.5.",l,
                 value = TRUE)[1]
    if(is.na(line)) return(NA)
    if(percent){
        as.numeric(gsub("[()]|[)]|%","",strsplit(line," ")[[1]][3]))
    } else { 
        as.integer(trimws(gsub("There are|,","",strsplit(line," SNPs")[[1]][1])))
    } 
}