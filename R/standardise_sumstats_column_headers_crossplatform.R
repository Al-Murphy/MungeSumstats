
#' Standardise the column headers in the Summary Statistics files (CROSSPLATFORM)
#'
#' Use a reference data table of common column header names (stored in sumstatsColHeaders.rda) convert them to a standard set, i.e. chromosome --> CHR
#'
#' This function does not check that all the required column headers are present
#'
#' The amended header is written directly back into the file
#'
#' @param first_line String containing the first line of the sumstats file
#'
#' @return The amended column headers (also the column headers will be written directly into the summary statistics file)
#'
#' @examples
#' #col_headers = standardise.sumstats.column.headers.crossplatform("~/Downloads/202040.assoc.tsv")
standardise_sumstats_column_headers_crossplatform <- function(first_line) {
  column_headers <- strsplit(first_line, "\t")[[1]]
  message("First line of summary statistics file: ")
  message(paste(first_line,collapse = "\t"))
  column_headers <- toupper(column_headers)

  for (headerI in seq_len(nrow(MungeSumstats::sumstatsColHeaders))) {
    un <- MungeSumstats::sumstatsColHeaders[headerI, 1]
    cr <- MungeSumstats::sumstatsColHeaders[headerI, 2]
    if (un %in% column_headers & (!cr %in% column_headers)) {
      column_headers <- gsub(sprintf("^%s$", un), cr, column_headers)
    }
  }
  new_first_line <- paste(column_headers, collapse = "\t")
  return(new_first_line)
}
