#' Ensure that no columns are duplicated
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
check_dup_col <- function(sumstats_file, path){
  ..notDup = NULL
  col_headers = names(sumstats_file)
  if(sum(duplicated(col_headers))>0){
    msg <- paste0("There are ",sum(duplicated(col_headers)),
                  " duplicated columns which will be removed.")
    message(msg)
    notDup <- which(!duplicated(col_headers))
    sumstats_file <- sumstats_file[, ..notDup]
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
