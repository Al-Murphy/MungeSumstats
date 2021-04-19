#' Ensure that there is at least one signed column in summary statistics file
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @return null
check_signed_col <- function(sumstats_file){
  col_headers <- names(sumstats_file)
  signed_stat_column_names <- c("Z","OR","BETA","LOG_ODDS","SIGNED_SUMSTAT")
  stp_msg <- paste0("ERROR: cannot find a column name representing signed ",
                    "statistic in GWAS sumstats file:\n",
                    "'Z','OR','BETA','LOG_ODDS','SIGNED_SUMSTAT'")
  if(sum(signed_stat_column_names %in% col_headers)<1){
    stop(stp_msg)
  }
}
