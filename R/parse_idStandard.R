#' Standardised IEU MRC OpenGWAS ID
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Character
#' @keywords internal
parse_idStandard <- function(l){
  line <- grep("Parsing .*. data column",l,
               value = TRUE)[1]
  if(is.na(line)) return(NA)
  trimws(strsplit(line," ")[[1]][2], whitespace = "'|[ ]")
}