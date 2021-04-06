#' Ensure that there is at least one signed column in summary statistics file
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @return null
check_signed_col <- function(sumstats_file){
  rows_of_data <- c(sumstats_file[1], sumstats_file[2])
  col_headers <- strsplit(rows_of_data[1], "\t")[[1]]
  print(paste0("Checking that there is at least one signed sumstats column ",
               "(eg.: Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT)"))
  signed_stat_column_names <- c("Z","OR","BETA","LOG_ODDS","SIGNED_SUMSTAT")
  if(sum(signed_stat_column_names %in% col_headers)<1 %in% col_headers){
    print("Header of file:")
    print(rows_of_data)
    stop(paste0("ERROR: cannot find a column name representing signed statisti",
                  "c in GWAS sumstats file. I.e. Z, OR, BETA"))
  }
}
