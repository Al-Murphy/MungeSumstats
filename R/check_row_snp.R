#' Ensure all rows have SNPs beginning with rs or SNP, drop those that don't
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
check_row_snp <- function(sumstats_file, path){
  #TODO use data table instead
  # All rows should start with either SNP or rs... if they don't drop them
  print("Dropping all rows that don't start with 'rs'")
  sumstats_file <- c(sumstats_file[1],sumstats_file[grepl("^rs",sumstats_file)])
  writeLines(text=sumstats_file, con = path)
  
  return(sumstats_file)
}  