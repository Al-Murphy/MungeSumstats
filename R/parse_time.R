#' Parse the total time taken the munge the file
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Character
#' @keywords internal
parse_time <- function(l){
    line <- grep("Done munging in",l,
                 value = TRUE)[1]
    if(is.na(line)) return(NA)
    as.numeric(trimws(strsplit(line," ")[[1]][4], whitespace = "'|[ ]"))
}