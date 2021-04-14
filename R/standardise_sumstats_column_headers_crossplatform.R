#' Standardise the column headers in the Summary Statistics files (CROSSPLATFORM)
#'
#' Use a reference data table of common column header names (stored in
#' sumstatsColHeaders) convert them to a standard set, i.e. chromosome -> CHR
#' This function does not check that all the required column headers are present
#' The amended header is written directly back into the file
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @examples
#' #col_headers = standardise.sumstats.column.headers.crossplatform("~/Downloads/202040.assoc.tsv")
standardise_sumstats_column_headers_crossplatform <-
  function(sumstats_file, path) {
  column_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  pretty_col_headers <- paste(sumstats_file[1],collapse = "\t")
  message("First line of summary statistics file: ")
  message(pretty_col_headers)
  column_headers <- toupper(column_headers)
  #Go through each and get correct spelling
  for (headerI in seq_len(nrow(MungeSumstats::sumstatsColHeaders))) {
    un <- MungeSumstats::sumstatsColHeaders[headerI, 1]
    cr <- MungeSumstats::sumstatsColHeaders[headerI, 2]
    if (un %in% column_headers & (!cr %in% column_headers)) {
      column_headers <- gsub(sprintf("^%s$", un), cr, column_headers)
    }
  }
  new_first_line <- paste(column_headers, collapse = "\t")
  sumstats_file[1] <- new_first_line
  writeLines(sumstats_file, con = path)

  return(sumstats_file)
}
