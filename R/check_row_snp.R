#' Ensure all rows have SNPs beginning with rs or SNP, drop those that don't
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @return The modified sumstats_file
#' @importFrom data.table fread
#' @importFrom data.table fwrite
check_row_snp <- function(sumstats_file, path){
  SNP = NULL
  # All rows should start with either SNP or rs... if they don't drop them
  #use data table for speed
  sumstats_dt <- data.table::fread(path)
  num_bad_ids <- nrow(sumstats_dt[!grep("^rs",SNP),])
  if(num_bad_ids>0){
    message(paste0(num_bad_ids, " SNPs",
                    " don't start with 'rs' and will be removed"))
    sumstats_dt <- sumstats_dt[grep("^rs",SNP),]
    #write new data
    data.table::fwrite(sumstats_dt,file=path,sep="\t")
    sumstats_file <- readLines(path)

    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
