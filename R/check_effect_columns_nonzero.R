#' Ensure that the standard error (se) is positive for all SNPs
#'
#' @inheritParams format_sumstats
#' @return list containing sumstats_dt, the modified summary statistics data 
#' table object
#' @keywords internal
check_effect_columns_nonzero <- 
    function(sumstats_dt, path, effect_columns_nonzero){
    col_headers <- names(sumstats_dt)
    effect_columns <- c("BETA","OR","LOG_ODDS","SIGNED_SUMSTAT")
    if(sum(effect_columns %in% col_headers)>=1 && effect_columns_nonzero){
        message("Filtering effect columns, ensuring none equal 0.")
        #filter to effect columns in the data
        effect_columns_dat <- effect_columns[effect_columns %in% col_headers]
        #ensure numeric
        sumstats_dt[,(effect_columns_dat):= lapply(.SD, as.numeric),
                        .SDcols = effect_columns_dat]
        #check if any equal 0 - use data table for speed
        bad_ids <-sumstats_dt[, Reduce(`|`, lapply(.SD, `==`, 0)),
                                .SDcols = effect_columns_dat]
        num_bad_ids <- sum(bad_ids)
        if(num_bad_ids>0){
            msg <- paste0(formatC(num_bad_ids,big.mark = ","), " SNPs",
                          " have effect values = 0 and will be removed")
            message(msg)
            sumstats_dt <- sumstats_dt[!bad_ids]
        }
        return(list("sumstats_dt"=sumstats_dt))
    }
    else{
        return(list("sumstats_dt"=sumstats_dt))
    }
}