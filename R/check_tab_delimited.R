#' Ensure that tabs separate rows
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @return The modified sumstats_file
check_tab_delimited <-function(sumstats_file){
  row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
  #check if readLines picked up headers as single char
  if (length(row_of_data) == 1){
    #check presence of space - space delimited
    if (length(grep(" ", row_of_data))!=0){
      print(paste0("WARNING: This GWAS sumstat file has space field separators",
                   " instead of tabs (unusual, not proper input for MAGMA). ",
                   "Temp file with corrected FS created and used instead."))
      sumstats_file <-
        gsub(pattern = " ", replacement = "\t", x = sumstats_file)
      return(sumstats_file)
    }
    #check presence of comma - comma delimited
    if (length(grep(",", row_of_data))!=0){
      print(paste0("WARNING: This GWAS sumstat file has comma field separators",
                   " instead of tabs (unusual, not proper input for MAGMA). ",
                   "Temp file with corrected FS created and used instead."))
      sumstats_file <-
        gsub(pattern = ",", replacement = "\t", x = sumstats_file)
      return(sumstats_file)
    }
  }
  else{
    return(sumstats_file)
  }
}
