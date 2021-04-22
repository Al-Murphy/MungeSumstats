#' Ensure all rows have SNPs beginning with rs or SNP, drop those that don't
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
#' @importFrom data.table fread
#' @importFrom data.table fwrite
check_row_snp <- function(sumstats_dt, path){
  SNP = NULL
  # All rows should start with either SNP or rs... if they don't drop them
  #use data table for speed
  num_bad_ids <- nrow(sumstats_dt[!grep("^rs",SNP),])
  if(num_bad_ids>0){
    msg <- paste0(num_bad_ids, " SNPs",
                  " don't start with 'rs' and will be removed")
    message(msg)
    sumstats_dt <- sumstats_dt[grep("^rs",SNP),]

    return(list("sumstats_dt"=sumstats_dt))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
