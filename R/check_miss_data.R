#' Ensure that only one model in GWAS sumstats or only one trait tested
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom stats complete.cases
check_miss_data <- function(sumstats_dt, path){
  #use data table for speed
  #check for rows missing data to be excluded
  if(nrow(sumstats_dt[!complete.cases(sumstats_dt),])>0){
    msg <- paste0("WARNING: ",
                  nrow(sumstats_dt[!complete.cases(sumstats_dt),]),
                  " rows in sumstats file are missing data and will ",
                  "be removed.")
    message(msg)
    sumstats_dt <- sumstats_dt[complete.cases(sumstats_dt)]
    
    return(list("sumstats_dt"=sumstats_dt))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
