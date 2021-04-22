#' Ensure that there is at least one signed column in summary statistics file
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @return null
#' @keywords internal
check_signed_col <- function(sumstats_dt){
  col_headers <- names(sumstats_dt)
  signed_stat_column_names <- c("Z","OR","BETA","LOG_ODDS","SIGNED_SUMSTAT")
  stp_msg <- paste0("ERROR: cannot find a column name representing signed ",
                    "statistic in GWAS sumstats file:\n",
                    "'Z','OR','BETA','LOG_ODDS','SIGNED_SUMSTAT'")
  if(sum(signed_stat_column_names %in% col_headers)<1){
    stop(stp_msg)
  }
}
