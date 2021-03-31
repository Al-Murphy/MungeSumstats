#' Ensure that no columns are duplicated
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted 
#' @return The modified sumstats_file
#' @importFrom utils read.table
#' @importFrom utils write.table
check_dup_col <- function(sumstats_file, path){
  print("Removing duplicated columns (if any).")
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]) 
  col_headers = strsplit(rows_of_data[1], "\t")[[1]]
  if(sum(duplicated(col_headers))>0){
    notDup <- which(!duplicated(col_headers))
    write.table(x=utils::read.table(path)[,notDup], file=path, sep="\t", 
                quote=FALSE, row.names = FALSE, col.names = FALSE) 
    sumstats_file <- readLines(path)
  }
  else{
    return(sumstats_file)
  }
}