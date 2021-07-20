#' Ensure that the standard error (se) is positive for all SNPs
#'
#' @inheritParams format_sumstats
#' @return list containing sumstats_dt, the modified summary statistics data 
#' table object
#' @keywords internal
check_pos_se <- function(sumstats_dt, path, pos_se){
    SE = NULL
    col_headers <- names(sumstats_dt)
    if("SE" %in% col_headers && pos_se){
        message("Filtering SNPs, ensuring Standard Error (SE) >0.")
        #use data table for speed
        num_bad_ids <- nrow(sumstats_dt[SE<=0,])
        if(num_bad_ids>0){
            msg <- paste0(formatC(num_bad_ids,big.mark = ","), " SNPs",
                          " have SE values <= 0 and will be removed")
            message(msg)
            sumstats_dt <- sumstats_dt[SE>0,]
        }
        return(list("sumstats_dt"=sumstats_dt))
    }
    else{
        return(list("sumstats_dt"=sumstats_dt))
    }
}