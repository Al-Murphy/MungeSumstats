#' Ensure that there is at least one signed column in summary statistics file
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @return null
check_signed_col <- function(sumstats_file){
  rows_of_data <- c(sumstats_file[1], sumstats_file[2])
  col_headers <- strsplit(rows_of_data[1], "\t")[[1]]
  signed_stat_column_names <- c("Z","OR","BETA","LOG_ODDS","SIGNED_SUMSTAT")
  stp_msg <- paste0("ERROR: cannot find a column name representing signed ",
                    "statistic in GWAS sumstats file:\n",
                    "'Z','OR','BETA','LOG_ODDS','SIGNED_SUMSTAT'")
  if(sum(signed_stat_column_names %in% col_headers)<1){
    stop(stp_msg)
  }
}
