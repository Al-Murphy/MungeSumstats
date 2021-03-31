#' Ensure all rows have SNPs beginning with rs or SNP, drop those that don't
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @importFrom data.table fread
#' @importFrom data.table fwrite
check_dup_snp <- function(sumstats_file, path){
  #TODO improve this
  # Try to remove duplicated RSIDs
  print("Removing duplicated RSIDs.")
  sumstats <- fread(path)
  if(sum(duplicated(sumstats[,1]))>0){
    notDup <- which(!duplicated(sumstats[,1]))
    notDupLines <- sumstats[notDup,]
    fwrite(notDupLines, file=path, sep="\t") 
    rm(notDupLines)
    gc()
  }
}