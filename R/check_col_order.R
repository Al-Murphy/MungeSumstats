#' Ensure that the first three columns are SNP, CHR, BP in that order
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table setcolorder
check_col_order <- function(sumstats_file, path){
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]) 
  col_headers <- strsplit(rows_of_data[1], "\t")[[1]]
  print(paste0("Checking that the first three column headers are SNP, CHR and",
                " BP in this order."))
  #Use data tables for speed
  if(!sum(col_headers[seq_len(3)]==c("SNP","CHR","BP"))==3){
    whichSNP <- which(col_headers=="SNP")[1]
    whichCHR <- which(col_headers=="CHR")[1]
    whichBP <- which(col_headers=="BP")[1]
    otherCols <- 
      setdiff(seq_len(length(col_headers)),c(whichSNP,whichCHR,whichBP))
    dt_sumstats <- data.table::fread(path)
    data.table::setcolorder(dt_sumstats, c(whichSNP,whichCHR,whichBP,otherCols))
    data.table::fwrite(x=dt_sumstats, file=path, sep="\t")
    sumstats_file <- readLines(path)
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}  