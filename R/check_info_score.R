#' Ensure all SNPs have info score above threshold
#'
#' @inheritParams format_sumstats
#' @return list containing sumstats_dt, the modified summary statistics data 
#' table object
#' @keywords internal
check_info_score <- function(sumstats_dt, path, INFO_filter){
  INFO = NULL
  col_headers <- names(sumstats_dt)
  if("INFO" %in% col_headers && INFO_filter>0){
    message("Filtering SNPs based on INFO score.")
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
