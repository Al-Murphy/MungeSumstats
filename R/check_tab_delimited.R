#' Ensure that tabs separate rows
#'
#' @param header The summary statistics file for the GWAS
#' @return The modified header
#' @keywords internal
check_tab_delimited <-function(header){
  row_of_data <- strsplit(header, "\t")[[1]]
  #check if readLines picked up headers as single char
  if (length(row_of_data) == 1){
    #check presence of space - space delimited
    if (length(grep(" ", row_of_data))!=0){
      msg <- paste0("WARNING: GWAS sumstat file has space field separators",
                    " instead of tabs (unusual, not proper input for MAGMA).\n",
                    "Temp file with corrected FS created and used instead.")
      message(msg)
      header <-
        gsub(pattern = " ", replacement = "\t", x = header)
      return(header)
    }
    #check presence of comma - comma delimited
    if (length(grep(",", row_of_data))!=0){
      msg <- paste0("WARNING: GWAS sumstat file has comma field separators",
                    " instead of tabs (unusual, not proper input for MAGMA).\n",
                      "Temp file with corrected FS created and used instead.")
      message(msg)
      header <-
        gsub(pattern = ",", replacement = "\t", x = header)
      return(header)
    }
  }
  else{
    return(header)
  }
}
