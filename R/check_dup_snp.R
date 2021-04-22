#' Ensure all rows have unique SNP IDs, drop those that don't
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table key
#' @importFrom data.table setkey
check_dup_snp <- function(sumstats_dt, path){
  SNP = NULL
  col_headers <- names(sumstats_dt)
  if("SNP" %in% col_headers){
    # Try to remove duplicated RSIDs
    data.table::setkey(sumstats_dt,SNP)
    dups <- duplicated(sumstats_dt, by = data.table::key(sumstats_dt))
    if(sum(dups)>0){
      dup_snps <- sumstats_dt$SNP[dups]
      msg <- paste0(paste(dup_snps, collapse = ", ")," RS IDs are duplicated ",
                    "in the sumstats file. These duplicates will be removed")
      message(msg)
      sumstats_dt <- unique(sumstats_dt, by = data.table::key(sumstats_dt))
  
      return(list("sumstats_dt"=sumstats_dt))
    }
  }
  return(list("sumstats_dt"=sumstats_dt))
}
