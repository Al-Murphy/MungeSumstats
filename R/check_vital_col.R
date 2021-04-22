#' Ensure that all necessary columns are in the summary statistics file
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @return null
#' @keywords internal
check_vital_col <- function(sumstats_dt){
  col_headers <- names(sumstats_dt)
  err_msg <-
    "Cannot find a %s column in GWAS sumstats. \nUse code such as '%s' to fix"
  for(key_column in c("SNP","CHR","BP","P","A1","A2")){
    code_example <- "sed -i '' '1s/p_value/P/' IQ.Sniekers.2017.txt"
    if(!key_column %in% col_headers){
      stop(sprintf(err_msg,key_column,code_example))
    }
  }
}
