#' Ensure that the N column are all integers
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table :=
check_n_int <- function(sumstats_file, path){
  N = N_tmp = NULL
  # Sometimes the N column is not all integers... so round it up
  col_headers <- strsplit(sumstats_file[1], "\t")[[1]]
  if("N" %in% col_headers) {
    sumstats_dt <- data.table::fread(path)
    if(is.integer(sumstats_dt$N)){ # check if any decimal places
      msg <- paste0("The sumstats N column is not all integers. Do you want ",
                      "MAGMA.celltyping to (if such instances exist) round ",
                      "them up? \n0 for NO, 1 for YES: ")
      choice <- 2
      while(!choice %in% c(0,1)){
        choice <- readline(msg)
        if(!choice %in% c(0,1)){
          message(paste0(choice," is not a valid option. Please try again"))
        }
      }
      if (as.logical(as.numeric(readline(msg)))) { #if user wants to correct
        sumstats_dt[,N:=round(N,0)]
        data.table::fwrite(x=sumstats_dt, file=path, sep="\t")
        sumstats_file <- readLines(path)

        return(sumstats_file)
      }
    }
  }
  return(sumstats_file)
}
