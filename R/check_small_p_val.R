#' Ensure that the p values are not 3e-400 or lower, if so set to 0
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table :=
check_small_p_val <- function(sumstats_file, path){
  P = NULL
  msg <- paste0("MAGMA cannot handle P-values as low as 3e-400. Do you want ",
                  "MAGMA.celltyping to \nconvert any (if) existing ones to ",
                  "zeroes? Type 0 for NO, 1 for YES: ")
  #Use data tables for speed
  sumstats <- fread(path)
  # MAGMA cannot handle P-values as low as 3e-400... so convert them to zeros
  if (as.logical(as.numeric(readline(msg)))) {
    rows_of_data <- c(sumstats_file[1], sumstats_file[2]) 
    col_headers <- strsplit(rows_of_data[1], "\t")[[1]]
    #TODO Note, I've not tested this since changing it from the original code... which was using gawk/sed
    sumstats[,P:=as.numeric(as.character(P))]
    fwrite(x=sumstats, file=path, sep="\t")
    sumstats_file <- readLines(path)
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}  