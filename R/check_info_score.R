#' Ensure all SNPs have info score above threshold
#'
#' @param sumstats_file The summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param INFO_filter integer The minimum value permissible of the imputation information score (if present in sumstatsfile). Default 0.9
#' @return The modified sumstats_file
check_info_score <- function(sumstats_file, path, INFO_filter){
  INFO = NULL
  col_headers <- names(sumstats_file)
  if("INFO" %in% col_headers && INFO_filter>0){
    #use data table for speed
    num_bad_ids <- nrow(sumstats_file[INFO<INFO_filter,])
    if(num_bad_ids>0){
      msg <- paste0(num_bad_ids, " SNPs",
                    " are below the INFO threshold of ",INFO_filter,
                    " and will be removed")
      message(msg)
      sumstats_file <- sumstats_file[INFO>=INFO_filter,]
    }
    return(sumstats_file)
  }
  else{
    return(sumstats_file)
  }
}
