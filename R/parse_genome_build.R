#' Genome build inferred from the summary statistics
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Character
#' @keywords internal
parse_genome_build <- function(l){
  line <- grep("Inferred genome build:",l,value = TRUE)[1]
  if(is.na(line)) return(NA)
  trimws(strsplit(line,":")[[1]][-1])
}