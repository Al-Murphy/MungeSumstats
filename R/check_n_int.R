#' Ensure that the N column are all integers
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param convert_n_int Binary, if N (the number of samples) is not an integer, should this be rounded? Default is TRUE.
#' @return The modified sumstats_file
#' @importFrom data.table :=
check_n_int <- function(sumstats_file, path, convert_n_int){
  N = N_tmp = NULL
  # Sometimes the N column is not all integers... so round it up
  col_headers <- names(sumstats_file)
  if("N" %in% col_headers) {
    if(!is.integer(sumstats_file$N)){ # check if any decimal places
      msg <- paste0("The sumstats N column is not all integers, this could ",
                      "effect downstream analysis.")
      if (convert_n_int) { #if user wants to correct
        msg2 <- paste0(msg,"These will be converted to integers.")
        message(msg2)
        sumstats_file[,N:=round(N,0)]

        return(sumstats_file)
      }
      else{
        msg2 <- paste0(msg,"These will NOT be converted to integers.")
        message(msg2)
      }
    }
  }
  return(sumstats_file)
}
