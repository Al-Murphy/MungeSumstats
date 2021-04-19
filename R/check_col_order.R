#' Ensure that the first three columns are SNP, CHR, BP in that order
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table setcolorder
check_col_order <- function(sumstats_file, path){
  col_headers <- names(sumstats_file)
  #Use data tables for speed
  if(!sum(col_headers[seq_len(3)]==c("SNP","CHR","BP"))==3){
    msg <- paste0("Reordering so first three column headers are SNP, CHR and",
                  " BP in this order.")
    message(msg)
    whichSNP <- which(col_headers=="SNP")[1]
    whichCHR <- which(col_headers=="CHR")[1]
    whichBP <- which(col_headers=="BP")[1]
    otherCols <-
      setdiff(seq_len(length(col_headers)),c(whichSNP,whichCHR,whichBP))
    data.table::setcolorder(sumstats_file, c(whichSNP,whichCHR,whichBP,otherCols))
    
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
