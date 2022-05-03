#' Parse number of SNPs dropped due to being on chrom X, Y or MT
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_dropped_chrom <- function(l){
    line <- grep("are on chromosomes X, Y, MT and will be removed",l,
                 value = TRUE)[1]
    if(is.na(line)) return(NA)
    as.integer(trimws(gsub(",","",strsplit(line," ")[[1]][1])))
}