#' Formatted example
#' 
#' Returns an example of summary stats that have had their column names 
#' already standardised with 
#' \link[MungeSumstats]{standardise_sumstats_column_headers_crossplatform}.
#' 
#' @param path Path to raw example file. Default to built-in dataset.
#' @param formatted Whether the column names should be formatted 
#' (default:\code{TRUE}).
#' @param sorted Whether the rows should be sorted by genomic coordinates
#' (default:\code{TRUE}).
#' @return \code{sumstats_dt}
#' @export
#' @examples 
#' sumstats_dt <- MungeSumstats::formatted_example()
formatted_example <- function(path=system.file("extdata", "eduAttainOkbay.txt",
                                               package = "MungeSumstats"),
                              formatted=TRUE,
                              sorted=TRUE){
    sumstats_dt <- suppressMessages(
        read_sumstats(path = path)
    )
    if(formatted){
        sumstats_dt <-
            standardise_sumstats_column_headers_crossplatform(
                sumstats_dt = sumstats_dt)$sumstats_dt 
    } else {
        if(sorted) {
            messager("Setting sorted=FALSE (required when formatted=FALSE).")
            sorted <- FALSE
        }
    }
    if(sorted){
        sumstats_dt <- sort_coords(sumstats_dt = sumstats_dt)
    }
    return(sumstats_dt)
}
