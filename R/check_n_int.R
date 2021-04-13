#' Ensure that the N column are all integers
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param convert_n_int Binary, if N (the number of samples) is not an integer, should this be rounded? Default is TRUE.
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table :=
check_n_int <- function(sumstats_file, path, convert_n_int){
  N = N_tmp = NULL
  # Sometimes the N column is not all integers... so round it up
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  if("N" %in% col_headers) {
    sumstats_dt <- data.table::fread(path)
    if(!is.integer(sumstats_dt$N)){ # check if any decimal places
      msg <- paste0("The sumstats N column is not all integers, this could ",
                      "effect downstream analysis.")
      if (convert_n_int) { #if user wants to correct
        message(paste0(msg,"These will be converted to integers."))
        sumstats_dt[,N:=round(N,0)]
        data.table::fwrite(x=sumstats_dt, file=path, sep="\t")
        sumstats_file <- readLines(path)

        return(sumstats_file)
      }
      else{
        message(paste0(msg,"These will NOT be converted to integers."))
      }
    }
  }
  return(sumstats_file)
}
