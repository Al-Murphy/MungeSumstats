#' Ensure all rows have unique SNP IDs, drop those that don't
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return null
#' @importFrom data.table key
#' @importFrom data.table setkey
check_dup_snp <- function(sumstats_file, path){
  SNP = NULL
  col_headers <- names(sumstats_file)
  if("SNP" %in% col_headers){
    # Try to remove duplicated RSIDs
    data.table::setkey(sumstats_file,SNP)
    dups <- duplicated(sumstats_file, by = data.table::key(sumstats_file))
    if(sum(dups)>0){
      dup_snps <- sumstats_file$SNP[dups]
      msg <- paste0(paste(dup_snps, collapse = ", ")," RS IDs are duplicated ",
                    "in the sumstats file. These duplicates will be removed")
      message(msg)
      sumstats_file <- unique(sumstats_file, by = data.table::key(sumstats_file))
  
      return(sumstats_file)
    }
  }
  return(sumstats_file)
}
