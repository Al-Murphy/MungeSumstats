#' Standardise the column headers in the Summary Statistics files
#'
#' Use a reference data table of common column header names (stored in
#' sumstatsColHeaders or user inputted mapping file) to convert them to a 
#' standard set, i.e. chromosome -> CHR. This function does not check that all 
#' the required column headers are present. The amended header is written 
#' directly back into the file
#'
#' @inheritParams format_sumstats
#' @return list containing sumstats_dt, the modified summary statistics data 
#' table object
#' @keywords internal
#' @importFrom data.table setnames
standardise_sumstats_column_headers_crossplatform <-
  function(sumstats_dt, path, mapping_file) { 
  message("Standardising column headers.")  
  message("First line of summary statistics file: ")
  msg <- paste0(names(sumstats_dt),split="\t")
  message(msg)
  #first make all column headers upper case
  data.table::setnames(sumstats_dt,toupper(names(sumstats_dt)))
  column_headers <- names(sumstats_dt)
  #load synonym mapping - internal data no loading
  #Go through each and get correct spelling 
  #allow for differing cases column names
  colnames(mapping_file) <- toupper(colnames(mapping_file))
  for (headerI in seq_len(nrow(mapping_file))) {
    un <- mapping_file[headerI, "UNCORRECTED"]
    cr <- mapping_file[headerI, "CORRECTED"]
    if (un %in% column_headers & (!cr %in% column_headers)) {
      data.table::setnames(sumstats_dt,un,cr)
    }
  }
  return(list("sumstats_dt"=sumstats_dt))
}
