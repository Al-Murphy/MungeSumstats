#' Ensure that all necessary columns are in the summary statistics file
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @return null
check_vital_col <- function(sumstats_file){
  rows_of_data <- c(sumstats_file[1], sumstats_file[2])
  col_headers <- strsplit(rows_of_data[1], "\t")[[1]]
  err_msg <-
    "Cannot find a %s column in GWAS sumstats. \nUse code such as '%s' to fix"
  for(key_column in c("SNP","CHR","BP","P","A1","A2")){
    code_example <- "sed -i '' '1s/p_value/P/' IQ.Sniekers.2017.txt"
    if(!key_column %in% col_headers){
      stop(sprintf(err_msg,key_column,code_example))
    }
  }
}
