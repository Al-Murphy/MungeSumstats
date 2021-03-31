#' Ensure that the p values are not 3e-400 or lower, if so set to 0
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table :=
check_n_int <- function(sumstats_file, path){
  N = NULL
  # Sometimes the N column is not all integers... so round it up
  #TODO - I edited this on 3rd April 2020 so it uses data.table... but don't have a dataset to check that it still works
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]) 
  col_headers = strsplit(rows_of_data[1], "\t")[[1]]
  if("N" %in% col_headers) {
    msg <- paste0("Sometimes the N column is not all integers. Do you want ",
                    "MAGMA.celltyping to (if such instances exist) round them ",
                    "up? \n0 for NO, 1 for YES: ")
    #whichN = which(col_headers %in% "N")
    if (as.logical(as.numeric(readline(msg)))) {
      #rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]
      sumstats <- fread(path) 
      #TODO # Note, I've not tested this since changing it from the original code... which was using gawk/sed
      sumstats[,N:=round(as.numeric(as.character(N)))]
      fwrite(x=sumstats, file=path, sep="\t") 
      sumstats_file <- readLines(path)
      #for (i in seq_along(sumstats[,which(col_headers=="N")])) {
      #  if (sumstats[i,which(col_headers=="N")]=="N") {next} # To skip the header.
      #  sumstats[i,which(col_headers=="N")] <- round(as.numeric(as.character(sumstats[i,which(col_headers=="N")]))) # This converts anything under 3e-400 to zeros.
      #}
      #write.table(x=sumstats, file=path, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE); sumstats_file <- readLines(path)
      return(sumstats_file)
    }
  }
  return(sumstats_file)

}