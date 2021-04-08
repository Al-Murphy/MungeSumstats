#' Ensure that only one model in GWAS sumstats or only one trait tested
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom stats complete.cases
check_miss_data <- function(sumstats_file, path){
  #use data table for speed
  dt_sumstats <- data.table::fread(path)
  #check for rows missing data to be excluded
  if(nrow(dt_sumstats[!complete.cases(dt_sumstats),])>0){
    message(paste0("WARNING: ",
                    nrow(dt_sumstats[!complete.cases(dt_sumstats),]),
                    " rows in sumstats file are missing data and will ",
                    "be removed."))
    dt_sumstats <- dt_sumstats[complete.cases(dt_sumstats)]
    data.table::fwrite(x=dt_sumstats, file=path, sep="\t")
    sumstats_file <- readLines(path)
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
