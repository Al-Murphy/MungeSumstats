#' Ensure that only one model in GWAS sumstats or only one trait tested
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom stats complete.cases
check_miss_data <- function(sumstats_file, path){
  #use data table for speed
  #check for rows missing data to be excluded
  if(nrow(sumstats_file[!complete.cases(sumstats_file),])>0){
    msg <- paste0("WARNING: ",
                  nrow(sumstats_file[!complete.cases(sumstats_file),]),
                  " rows in sumstats file are missing data and will ",
                  "be removed.")
    message(msg)
    sumstats_file <- sumstats_file[complete.cases(sumstats_file)]
    
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
