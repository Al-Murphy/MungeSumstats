#' Standardise the column headers in the Summary Statistics files (CROSSPLATFORM)
#'
#' Use a reference data table of common column header names (stored in
#' sumstatsColHeaders) convert them to a standard set, i.e. chromosome -> CHR
#' This function does not check that all the required column headers are present
#' The amended header is written directly back into the file
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table setnames
standardise_sumstats_column_headers_crossplatform <-
  function(sumstats_dt, path) {
  message("First line of summary statistics file: ")
  msg <- paste0(names(sumstats_dt),split="\t")
  message(msg)
  #first make all column headers upper case
  data.table::setnames(sumstats_dt,toupper(names(sumstats_dt)))
  column_headers <- names(sumstats_dt)
  #load synonym mapping - internal data no loading
  #Go through each and get correct spelling
  for (headerI in seq_len(nrow(sumstatsColHeaders))) {
    un <- sumstatsColHeaders[headerI, 1]
    cr <- sumstatsColHeaders[headerI, 2]
    if (un %in% column_headers & (!cr %in% column_headers)) {
      data.table::setnames(sumstats_dt,un,cr)
    }
  }
  return(list("sumstats_dt"=sumstats_dt))
}
