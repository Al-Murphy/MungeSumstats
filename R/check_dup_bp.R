#' Ensure all rows have unique positions, drop those that don't
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table key
#' @importFrom data.table setkey
check_dup_bp <- function(sumstats_dt, path){
  BP = NULL
  col_headers <- names(sumstats_dt)
  if("BP" %in% col_headers){
    message("Checking for SNPs with duplicated base-pair positions")
    # Try to remove duplicated Positions
    data.table::setkey(sumstats_dt,BP)
    dups <- duplicated(sumstats_dt, by = data.table::key(sumstats_dt))
    if(sum(dups)>0){
      msg <- paste0(formatC(sum(dups),big.mark = ",")," base-pair positions are",
                    " duplicated in the sumstats file. These duplicates will",
                    " be removed")
      message(msg)
      sumstats_dt <- unique(sumstats_dt, by = data.table::key(sumstats_dt))
      
      return(list("sumstats_dt"=sumstats_dt))
    }
  }
  return(list("sumstats_dt"=sumstats_dt))
}
