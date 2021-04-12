#' Ensure all rows have SNPs beginning with rs or SNP, drop those that don't
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return null
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom data.table key
#' @importFrom data.table setkey
check_dup_snp <- function(sumstats_file, path){
  SNP = NULL
  # Try to remove duplicated RSIDs
  sumstats_dt <- data.table::fread(path)
  data.table::setkey(sumstats_dt,SNP)
  dups <- duplicated(sumstats_dt, by = data.table::key(sumstats_dt))
  if(sum(dups)>0){
    message(paste0(dups," RS IDs are duplicated ",
                    "in the sumstats file. These duplicates will be removed"))
    sumstats_dt <- unique(sumstats_dt, by = data.table::key(sumstats_dt))
    data.table::fwrite(sumstats_dt, file=path, sep="\t")
    sumstats_file <- readLines(path)

    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
