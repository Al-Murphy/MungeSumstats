#' Ensure all SNPs have info score above threshold
#'
#' @param sumstats_dt data table obj of the summary statistics file for the GWAS
#' @param path Filepath for the summary statistics file to be formatted
#' @param INFO_filter integer The minimum value permissible of the imputation information score (if present in sumstatsfile). Default 0.9
#' @return list containing sumstats_dt, the modified summary statistics data table object
#' @keywords internal
check_info_score <- function(sumstats_dt, path, INFO_filter){
  INFO = NULL
  col_headers <- names(sumstats_dt)
  if("INFO" %in% col_headers && INFO_filter>0){
    #use data table for speed
    num_bad_ids <- nrow(sumstats_dt[INFO<INFO_filter,])
    if(num_bad_ids>0){
      msg <- paste0(formatC(num_bad_ids,big.mark = ","), " SNPs",
                    " are below the INFO threshold of ",INFO_filter,
                    " and will be removed")
      message(msg)
      sumstats_dt <- sumstats_dt[INFO>=INFO_filter,]
    }
    return(list("sumstats_dt"=sumstats_dt))
  }
  else{
    return(list("sumstats_dt"=sumstats_dt))
  }
}
