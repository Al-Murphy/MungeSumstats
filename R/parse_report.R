#' Parse "Summary statistics report" metrics 
#' 
#' Support function for \link[MungeSumstats]{parse_logs}.
#' 
#' @param l Lines of text from log file.
#' 
#' @return Numeric
#' @keywords internal
parse_report <- function(l, 
                         entry = 1,
                         line = 1){
    report_lines <- grep("Summary statistics report:",l)
    if(length(report_lines)<2 && entry!=1) return(NA)
    as.integer(
        gsub(",","",
             strsplit(
                 trimws(
                     l[report_lines+line][entry]
                 ),
                 " ")[[1]][2]
        )
    )
}